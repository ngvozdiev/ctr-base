#include <gflags/gflags.h>
#include <stddef.h>
#include <SFML/Graphics.hpp>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <deque>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/animator.h"
#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/htsim/packet.h"
#include "ncode_common/src/htsim/queue.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/strutil.h"
#include "../common.h"
#include "../net_mock.h"
#include "../pcap_data.h"
#include "manual_network.h"
#include "sfml.h"

DEFINE_string(pcap_trace_store, "trace_store.pb",
              "A file with information about .pcap traces");

// File to use for all fonts.
static constexpr char kDefaultFontFile[] = "DejaVuSans.ttf";

static void BootstrapNetwork(const ctr::PcapTraceStore& trace_store,
                             const std::vector<ctr::AggregateId>& aggregates,
                             ctr::manual::NetworkContainer* network_container,
                             ctr::MockSimDeviceFactory* device_factory) {
  const nc::net::GraphStorage* graph = network_container->graph();

  // Bin sequences will be allocated at random.
  std::map<ctr::AggregateId, std::unique_ptr<ctr::BinSequence>> sequences;
  std::vector<const ctr::PcapDataTrace*> all_traces = trace_store.AllTraces();
  CHECK(all_traces.size() >= aggregates.size());
  for (size_t i = 0; i < aggregates.size(); ++i) {
    const ctr::AggregateId& id = aggregates[i];
    const ctr::PcapDataTrace* trace = all_traces[i];
    std::vector<ctr::BinSequence::TraceAndSlice> slices =
        trace->TracesAndSlices(trace->AllSlices());

    auto bin_sequence = nc::make_unique<ctr::BinSequence>(slices);
    nc::htsim::MatchRuleKey key_for_aggregate =
        network_container->AddAggregate(id);

    // Will add the bin sequence to the source device.
    const std::string& src_device_id = graph->GetNode(id.src())->id();
    device_factory->AddBinSequence(src_device_id, key_for_aggregate,
                                   std::move(bin_sequence));
  }

  network_container->AddElementsFromGraph(device_factory, nullptr, nullptr);
  device_factory->Init();
}

class QueueMonitor : public nc::EventConsumer {
 public:
  QueueMonitor(const nc::htsim::Queue* queue, ctr::sfml::VisualLink* visual,
               nc::EventQueue* event_queue)
      : nc::EventConsumer(nc::StrCat("QueueMonitor", queue->id()), event_queue),
        period_(visual->time_step()),
        queue_(queue),
        visual_(visual) {
    EnqueueIn(event_queue->ToTime(period_));
  }

  void HandleEvent() override {
    Record();
    EnqueueIn(event_queue()->ToTime(period_));
  }

 private:
  void Record() {
    const nc::htsim::QueueStats& stats = queue_->GetStats();
    if (stats.tag_to_bytes_seen.empty()) {
      return;
    }

    std::map<nc::htsim::PacketTag, double> to_update;
    for (const auto& tag_and_bytes_seen : stats.tag_to_bytes_seen) {
      nc::htsim::PacketTag tag = tag_and_bytes_seen.first;
      uint64_t bytes_seen = tag_and_bytes_seen.second;
      uint64_t delta_bytes = bytes_seen - prev_bytes_seen_[tag];

      // This is bytes per period. Need to convert to Mbps.
      double periods_per_sec = 1000.0 / period_.count();
      double bits_per_sec = delta_bytes * 8 * periods_per_sec;
      double delta_Mbps = bits_per_sec / 1000000.0;
      to_update.emplace(tag, delta_Mbps);
    }

    visual_->AddValues(to_update);
    prev_bytes_seen_ = stats.tag_to_bytes_seen;
  }

  std::chrono::milliseconds period_;
  std::map<nc::htsim::PacketTag, uint64_t> prev_bytes_seen_;
  const nc::htsim::Queue* queue_;
  ctr::sfml::VisualLink* visual_;
};

std::tuple<ctr::sfml::VisualLinkStyle, ctr::sfml::NodeStyle> GetStyle() {
  ctr::sfml::NodeStyle node_style;
  sf::Color node_fill_color = sf::Color(203, 225, 243);
  node_style.fill = node_fill_color;
  node_style.outline = sf::Color::Black;
  node_style.outline_thickness = 2;

  ctr::sfml::VisualLinkStyle link_style;
  sf::Color link_background_color = sf::Color(230, 230, 250);
  link_style.background_color = link_background_color;

  ctr::sfml::LineStyle base_line_style;
  base_line_style.color = sf::Color::Black;
  base_line_style.thickness = 4;

  ctr::sfml::LineStyle top_dash_style;
  top_dash_style.color = sf::Color::Black;
  top_dash_style.thickness = 2;

  ctr::sfml::DashedLineStyle top_line_style;
  top_line_style.dash_len = 10;
  top_line_style.skip_len = 5;
  top_line_style.line_style = top_dash_style;

  link_style.base_line_style = base_line_style;
  link_style.top_line_style = top_line_style;

  return std::make_tuple(link_style, node_style);
}

static constexpr char kFrankfurt[] = "Frankfurt";
static constexpr char kBudapest[] = "Budapest";
static constexpr char kVienna[] = "Vienna";
static constexpr char kInternet[] = "Internet";
static constexpr char kB2V[] = "B2V";

static nc::net::GraphBuilder GetGraph() {
  nc::net::GraphBuilder builder;
  builder.AddLink({kBudapest, kVienna,
                   nc::net::Bandwidth::FromGBitsPerSecond(1.8),
                   std::chrono::milliseconds(500)});
  builder.AddLink({kVienna, kBudapest,
                   nc::net::Bandwidth::FromGBitsPerSecond(1.8),
                   std::chrono::milliseconds(500)});

  builder.AddLink({kFrankfurt, kBudapest,
                   nc::net::Bandwidth::FromGBitsPerSecond(3.6),
                   std::chrono::seconds(1)});
  builder.AddLink({kBudapest, kFrankfurt,
                   nc::net::Bandwidth::FromGBitsPerSecond(3.6),
                   std::chrono::seconds(1)});

  builder.AddLink({kFrankfurt, kVienna,
                   nc::net::Bandwidth::FromGBitsPerSecond(3.6),
                   std::chrono::seconds(1)});
  builder.AddLink({kVienna, kFrankfurt,
                   nc::net::Bandwidth::FromGBitsPerSecond(3.6),
                   std::chrono::seconds(1)});

  builder.AddLink({kInternet, kFrankfurt,
                   nc::net::Bandwidth::FromGBitsPerSecond(10),
                   std::chrono::seconds(1)});
  builder.AddLink({kB2V, kBudapest, nc::net::Bandwidth::FromGBitsPerSecond(10),
                   std::chrono::seconds(1)});

  return builder;
}

using PathPtr = std::unique_ptr<nc::net::Walk>;

static std::tuple<PathPtr, PathPtr, PathPtr> SetUpRouting(
    const ctr::PcapTraceStore& trace_store,
    ctr::manual::NetworkContainer* container,
    ctr::MockSimDeviceFactory* device_factory,
    ctr::sfml::VisualLinkStyle* link_style) {
  const nc::net::GraphStorage* graph = container->graph();

  ctr::AggregateId budapest_to_vienna(graph->NodeFromStringOrDie(kB2V),
                                      graph->NodeFromStringOrDie(kVienna));
  ctr::AggregateId frankfurt_to_vienna(graph->NodeFromStringOrDie(kInternet),
                                       graph->NodeFromStringOrDie(kVienna));

  BootstrapNetwork(trace_store, {budapest_to_vienna, frankfurt_to_vienna},
                   container, device_factory);

  auto p1 = graph->WalkFromStringOrDie(
      "[B2V->Budapest, Budapest->Frankfurt, Frankfurt->Vienna]");
  auto p2 = graph->WalkFromStringOrDie("[B2V->Budapest, Budapest->Vienna]");
  container->InstallPaths(budapest_to_vienna,
                          {{p1.get(), 0.5}, {p2.get(), 0.5}});

  auto p3 =
      graph->WalkFromStringOrDie("[Internet->Frankfurt, Frankfurt->Vienna]");
  container->InstallPaths(frankfurt_to_vienna, {{p3.get(), 1.0}});

  nc::htsim::PacketTag p1tag = container->TagForPathOrDie(*p1);
  nc::htsim::PacketTag p2tag = container->TagForPathOrDie(*p2);
  nc::htsim::PacketTag p3tag = container->TagForPathOrDie(*p3);

  link_style->plot_colors[p1tag] = sf::Color(85, 107, 47);
  link_style->plot_colors[p2tag] = sf::Color(85, 107, 47);
  link_style->plot_colors[p3tag] = sf::Color(205, 92, 92);

  return std::make_tuple(std::move(p1), std::move(p2), std::move(p3));
}

void RebalanceRoutes(double fraction_via_frankfurt,
                     const nc::net::Walk* via_frankfurt,
                     const nc::net::Walk* direct,
                     ctr::manual::NetworkContainer* container) {
  const nc::net::GraphStorage* graph = container->graph();
  ctr::AggregateId budapest_to_vienna(graph->NodeFromStringOrDie(kB2V),
                                      graph->NodeFromStringOrDie(kVienna));
  container->InstallPaths(budapest_to_vienna,
                          {{via_frankfurt, fraction_via_frankfurt},
                           {direct, 1 - fraction_via_frankfurt}});
}

// A combination of time in seconds and a location.
using TimeAndLocation = std::pair<double, sf::Vector2f>;

// A combination of time in seconds and a zoom value.
using TimeAndZoom = std::pair<double, float>;

static void AnimateViewPosition(
    const std::vector<TimeAndLocation>& frames, sf::View* view,
    nc::htsim::AnimationContainer* animation_container) {
  std::vector<nc::htsim::KeyFrame> x_key_frames;
  std::vector<nc::htsim::KeyFrame> y_key_frames;

  // Will always start at the current view.
  x_key_frames.emplace_back(std::chrono::milliseconds(0), view->getCenter().x);
  y_key_frames.emplace_back(std::chrono::milliseconds(0), view->getCenter().y);
  for (const auto& frame : frames) {
    std::chrono::milliseconds at(static_cast<size_t>(frame.first * 1000));
    x_key_frames.emplace_back(at, frame.second.x);
    y_key_frames.emplace_back(at, frame.second.y);
  }

  auto view_x_animator = nc::make_unique<nc::htsim::LinearAnimator>(
      x_key_frames, false, [view](double value) {
        float y = view->getCenter().y;
        view->setCenter(value, y);
      });

  auto view_y_animator = nc::make_unique<nc::htsim::LinearAnimator>(
      y_key_frames, false, [view](double value) {
        float x = view->getCenter().x;
        view->setCenter(x, value);
      });

  animation_container->AddAnimator(std::move(view_x_animator));
  animation_container->AddAnimator(std::move(view_y_animator));
}

static void AnimateViewZoom(
    const std::vector<TimeAndZoom>& frames, sf::View* view,
    nc::htsim::AnimationContainer* animation_container) {
  std::vector<nc::htsim::KeyFrame> key_frames;

  // Will always start at the current view.
  key_frames.emplace_back(std::chrono::milliseconds(0), 1.0);
  for (const auto& frame : frames) {
    std::chrono::milliseconds at(static_cast<size_t>(frame.first * 1000));
    key_frames.emplace_back(at, frame.second);
  }

  auto view_animator = nc::make_unique<nc::htsim::LinearAnimator>(
      key_frames, false, [view](double value) { view->zoom(value); });
  animation_container->AddAnimator(std::move(view_animator));
}

int main() {
  sf::ContextSettings settings;
  settings.antialiasingLevel = 4;

  sf::RenderWindow window(sf::VideoMode(800, 600), "SFML works!",
                          sf::Style::Default, settings);
  window.clear(sf::Color::White);
  window.setFramerateLimit(60);

  ctr::sfml::VisualLinkStyle link_style;
  ctr::sfml::NodeStyle node_style;
  std::tie(link_style, node_style) = GetStyle();

  ctr::PcapTraceStore trace_store(FLAGS_pcap_trace_store);
  nc::net::GraphBuilder graph_builder = GetGraph();
  nc::net::GraphStorage graph(graph_builder);
  nc::SimTimeEventQueue event_queue;

  ctr::manual::NetworkContainer network_container(
      std::chrono::milliseconds(100), &graph, &event_queue);
  ctr::MockSimDeviceFactory device_factory(
      ctr::manual::NetworkContainer::kDefaultEnterPort, &event_queue);

  PathPtr p1, p2, p3;
  std::tie(p1, p2, p3) = SetUpRouting(trace_store, &network_container,
                                      &device_factory, &link_style);

  sf::Vector2f budapest_location(600, 300);
  sf::Vector2f vienna_location(500, 400);
  sf::Vector2f frankfurt_location(150, 150);

  const nc::net::GraphLinkBase* budapest_frankfurt_link;
  nc::htsim::Queue* budapest_frankfurt_queue;
  std::tie(budapest_frankfurt_link, budapest_frankfurt_queue) =
      network_container.GetLink(kBudapest, kFrankfurt);
  budapest_frankfurt_queue->set_tags_in_stats(true);
  ctr::sfml::VisualLink vlink_budapest_frankfurt(
      budapest_location, frankfurt_location, 50, true,
      std::chrono::milliseconds(2), *budapest_frankfurt_link, link_style);
  QueueMonitor monitor_one(budapest_frankfurt_queue, &vlink_budapest_frankfurt,
                           &event_queue);

  const nc::net::GraphLinkBase* budapest_vienna_link;
  nc::htsim::Queue* budapest_vienna_queue;
  std::tie(budapest_vienna_link, budapest_vienna_queue) =
      network_container.GetLink(kBudapest, kVienna);
  budapest_vienna_queue->set_tags_in_stats(true);
  ctr::sfml::VisualLink vlink_budapest_vienna(
      budapest_location, vienna_location, 25, false,
      std::chrono::milliseconds(2), *budapest_vienna_link, link_style);
  QueueMonitor monitor_two(budapest_vienna_queue, &vlink_budapest_vienna,
                           &event_queue);

  const nc::net::GraphLinkBase* frankfurt_vienna_link;
  nc::htsim::Queue* frankfurt_vienna_queue;
  std::tie(frankfurt_vienna_link, frankfurt_vienna_queue) =
      network_container.GetLink(kFrankfurt, kVienna);
  frankfurt_vienna_queue->set_tags_in_stats(true);
  ctr::sfml::VisualLink vlink_frankfurt_vienna(
      frankfurt_location, vienna_location, 50, true,
      std::chrono::milliseconds(2), *frankfurt_vienna_link, link_style);
  QueueMonitor monitor_three(frankfurt_vienna_queue, &vlink_frankfurt_vienna,
                             &event_queue);

  ctr::sfml::Node budapest_node(20, budapest_location, node_style);
  ctr::sfml::Node vienna_node(20, vienna_location, node_style);
  ctr::sfml::Node frankfurt_node(20, frankfurt_location, node_style);

  sf::View view = window.getDefaultView();

  // A separate event queue for view changes.
  nc::SimTimeEventQueue view_event_queue;
  nc::htsim::AnimationContainer animation_container(
      "AnimationContainer", std::chrono::milliseconds(10), &view_event_queue);
  AnimateViewPosition({{10.0, sf::Vector2f(300, 300)}}, &view,
                      &animation_container);

  double rt_factor = 1.0;
  double via_frankfurt = 0.5;
  double step = 0.1;
  sf::Clock clock;
  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) {
        window.close();
      }

      if (event.key.code == sf::Keyboard::Up) {
        via_frankfurt += step;
      } else if (event.key.code == sf::Keyboard::Down) {
        via_frankfurt -= step;
      } else if (event.key.code == sf::Keyboard::Left) {
        rt_factor -= 0.1;
      } else if (event.key.code == sf::Keyboard::Right) {
        rt_factor += 0.1;
      }

      via_frankfurt = std::max(0.001, via_frankfurt);
      via_frankfurt = std::min(1.0, via_frankfurt);

      rt_factor = std::max(0.001, rt_factor);
      rt_factor = std::min(1.0, rt_factor);
    }

    sf::Time time = clock.restart();
    size_t elapsed = time.asMilliseconds();
    elapsed = std::min(elapsed, static_cast<size_t>(20));

    std::chrono::milliseconds delta_ms(elapsed);
    std::chrono::milliseconds delta_ms_after_factor(
        static_cast<size_t>(elapsed * rt_factor));

    RebalanceRoutes(via_frankfurt, p1.get(), p2.get(), &network_container);
    view_event_queue.RunAndStopIn(delta_ms);
    event_queue.RunAndStopIn(delta_ms_after_factor);

    window.setView(view);
    window.clear(sf::Color::White);
    window.draw(vlink_budapest_frankfurt);
    window.draw(vlink_budapest_vienna);
    window.draw(vlink_frankfurt_vienna);
    window.draw(budapest_node);
    window.draw(vienna_node);
    window.draw(frankfurt_node);
    window.display();
  }

  return 0;
}
