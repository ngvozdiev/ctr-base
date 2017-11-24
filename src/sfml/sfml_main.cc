#include <gflags/gflags.h>
#include <stddef.h>
#include <SFML/Graphics.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/strutil.h"
#include "../common.h"
#include "../controller.h"
#include "../net_mock.h"
#include "manual_network.h"

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
      : nc::EventConsumer("QueueMonitor", event_queue),
        period_(visual->time_step()),
        prev_bytes_seen_(0),
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
    uint64_t delta_bytes = stats.bytes_seen - prev_bytes_seen_;

    // This is bytes per period. Need to convert to Mbps.
    double periods_per_sec = 1000.0 / period_.count();
    double bits_per_sec = delta_bytes * 8 * periods_per_sec;
    double delta_Mbps = bits_per_sec / 1000000.0;
    visual_->AddValue(delta_Mbps);
    prev_bytes_seen_ = stats.bytes_seen;
  }

  std::chrono::milliseconds period_;
  uint64_t prev_bytes_seen_;
  const nc::htsim::Queue* queue_;
  ctr::sfml::VisualLink* visual_;
};

std::tuple<ctr::sfml::VisualLinkStyle, ctr::sfml::NodeStyle> GetStyle() {
  ctr::sfml::NodeStyle node_style;
  sf::Color node_fill_color = sf::Color(203, 225, 243, 1);
  node_style.fill = node_fill_color;
  node_style.outline = sf::Color::Black;
  node_style.outline_thickness = 2;

  ctr::sfml::VisualLinkStyle link_style;
  sf::Color link_background_color = sf::Color(221, 206, 255);
  sf::Color link_plot_color = sf::Color(150, 105, 255);
  link_style.background_color = link_background_color;
  link_style.plot_color = link_plot_color;

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

  return {link_style, node_style};
}

int main() {
  sf::ContextSettings settings;
  settings.antialiasingLevel = 4;

  sf::RenderWindow window(sf::VideoMode(800, 600), "SFML works!",
                          sf::Style::Default, settings);
  window.clear(sf::Color::White);
  window.setFramerateLimit(60);

  ctr::sfml::VisualLinkStyle link_style;
  ctr::sfml::VisualLinkStyle node_style;
  std::tie(link_style, node_style) = GetStyle();

  ctr::PcapTraceStore trace_store(FLAGS_pcap_trace_store);
  nc::net::GraphBuilder graph_builder = nc::net::GenerateFullGraph(
      2, nc::net::Bandwidth::FromGBitsPerSecond(5), std::chrono::seconds(10));
  nc::net::GraphStorage graph(graph_builder);
  nc::SimTimeEventQueue event_queue;

  ctr::manual::NetworkContainer network_container(
      std::chrono::milliseconds(100), &graph, &event_queue);
  ctr::MockSimDeviceFactory device_factory(
      ctr::manual::NetworkContainer::kDefaultEnterPort, &event_queue);

  // Will only have one aggregate.
  ctr::AggregateId aggregate_id(graph.NodeFromStringOrDie("N0"),
                                graph.NodeFromStringOrDie("N1"));
  BootstrapNetwork(trace_store, {aggregate_id}, &network_container,
                   &device_factory);

  const nc::net::GraphLinkBase* link;
  const nc::htsim::Queue* queue;
  std::tie(link, queue) = network_container.GetLink("N0", "N1");
  ctr::sfml::VisualLink vlink(sf::Vector2f(100, 100), sf::Vector2f(500, 500),
                              50, std::chrono::milliseconds(20), *link,
                              link_style);

  auto path = graph.WalkFromStringOrDie("[N0->N1]");
  network_container.InstallPaths(aggregate_id, {{path.get(), 1.0}});

  std::vector<const nc::htsim::Queue*> all_queues = network_container.queues();
  QueueMonitor monitor(all_queues[0], &vlink, &event_queue);

  ctr::sfml::Node n1(20, 100, 100, node_style);
  ctr::sfml::Node n2(20, 500, 500, node_style);

  sf::Clock clock;
  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) window.close();
    }

    sf::Time time = clock.restart();
    std::chrono::milliseconds delta_ms(time.asMilliseconds());
    event_queue.RunAndStopIn(delta_ms);

    window.clear(sf::Color::White);
    window.draw(vlink);
    window.draw(n1);
    window.draw(n2);
    window.display();
  }

  return 0;
}
