#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/strutil.h"
#include "common.h"
#include "controller.h"
#include "mean_est/mean_est.h"
#include "metrics/metrics.h"
#include "net_instrument.h"
#include "net_mock.h"
#include "opt/ctr.h"
#include "opt/opt.h"
#include "opt/path_provider.h"
#include "pcap_data.h"
#include "routing_system.h"
#include "tldr.h"

DEFINE_string(topology, "", "A file with a topology");
DEFINE_string(traffic_matrix, "", "A file with a traffic matrix");
DEFINE_string(pcap_trace_store, "trace_store.pb",
              "A file with information about .pcap traces");
DEFINE_uint64(period_duration_ms, 60000, "Length of the period");
DEFINE_uint64(history_bin_size_ms, 100, "How big each history bin is");
DEFINE_double(decay_factor, 0.0, "How quickly to decay prediction");
DEFINE_double(link_capacity_scale, 1.0, "By how much to scale all links");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_string(opt, "CTR", "The optimizer to use");

template <typename T>
static void PrintTimeDiff(std::ostream& out, T chrono_diff) {
  namespace sc = std::chrono;
  auto diff = sc::duration_cast<sc::milliseconds>(chrono_diff).count();
  auto const msecs = diff % 1000;
  diff /= 1000;
  auto const secs = diff % 60;
  diff /= 60;
  auto const mins = diff % 60;
  diff /= 60;
  auto const hours = diff % 24;
  diff /= 24;
  auto const days = diff;

  bool printed_earlier = false;
  if (days >= 1) {
    printed_earlier = true;
    out << days << (1 != days ? " days" : " day") << ' ';
  }
  if (printed_earlier || hours >= 1) {
    printed_earlier = true;
    out << hours << (1 != hours ? " hours" : " hour") << ' ';
  }
  if (printed_earlier || mins >= 1) {
    printed_earlier = true;
    out << mins << (1 != mins ? " minutes" : " minute") << ' ';
  }
  if (printed_earlier || secs >= 1) {
    printed_earlier = true;
    out << secs << (1 != secs ? " seconds" : " second") << ' ';
  }
  if (printed_earlier || msecs >= 1) {
    printed_earlier = true;
    out << msecs << (1 != msecs ? " milliseconds" : " millisecond");
  }
}

static std::chrono::milliseconds TimeNow() {
  return std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::high_resolution_clock::now().time_since_epoch());
}

class ProgressIndicator : public nc::EventConsumer {
 public:
  ProgressIndicator(std::chrono::milliseconds update_period,
                    nc::EventQueue* event_queue)
      : nc::EventConsumer("ProgressIndicator", event_queue),
        period_(event_queue->ToTime(update_period)),
        init_real_time_(TimeNow()) {
    EnqueueIn(period_);
  }

  void HandleEvent() override {
    double progress = event_queue()->Progress();
    CHECK(progress >= 0 && progress <= 1.0);
    std::cout << "\rProgress: " << std::setprecision(3) << (progress * 100.0)
              << "% ";

    auto real_time_delta = TimeNow() - init_real_time_;
    if (real_time_delta.count() > 0) {
      std::cout << "time remaining: ";
      auto remaining = std::chrono::milliseconds(static_cast<uint64_t>(
          real_time_delta.count() / progress * (1 - progress)));

      PrintTimeDiff(std::cout, remaining);
      std::cout << "                ";
    }

    std::cout << std::flush;
    EnqueueIn(period_);
  }

 private:
  // How often (in simulated time) to update the indicator.
  nc::EventQueueTime period_;

  // The initial time.
  std::chrono::milliseconds init_real_time_;
};

// A global variable that will keep a reference to the event queue, useful for
// logging, as the logging handler only accepts a C-style function pointer.
static nc::EventQueue* event_queue_global_ptr;

void CustomLogHandler(nc::LogLevel level, const char* filename, int line,
                      const std::string& message, nc::LogColor color) {
  uint64_t time_as_ms = event_queue_global_ptr->TimeToRawMillis(
      event_queue_global_ptr->CurrentTime());
  std::string new_message = nc::StrCat(time_as_ms, "ms ", message);
  nc::DefaultLogHandler(level, filename, line, new_message, color);
}

// Returns an aggregate history that will have a given mean rate.
static ctr::AggregateHistory GetDummyHistory(nc::net::Bandwidth rate,
                                             size_t flow_count) {
  size_t bin_size_ms = FLAGS_history_bin_size_ms;
  size_t init_window_ms = FLAGS_period_duration_ms;
  size_t bins_count = init_window_ms / bin_size_ms;

  double bins_in_second = 1000.0 / bin_size_ms;
  CHECK(bins_in_second > 0);
  double bytes_per_bin = (rate.bps() / 8.0) / bins_in_second;
  CHECK(bytes_per_bin > 0);

  std::vector<uint64_t> bins;
  for (size_t i = 0; i < bins_count; ++i) {
    bins.emplace_back(bytes_per_bin);
  }

  return {bins, std::chrono::milliseconds(bin_size_ms), flow_count};
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::metrics::InitMetrics();

  // The event queue.
  nc::SimTimeEventQueue event_queue;
  event_queue_global_ptr = &event_queue;
  nc::SetLogHandler(CustomLogHandler);

  auto timestamp_provider =
      ::nc::make_unique<nc::metrics::SimTimestampProvider>(&event_queue);
  nc::metrics::DefaultMetricManager()->set_timestamp_provider(
      std::move(timestamp_provider));

  CHECK(!FLAGS_topology.empty());
  CHECK(!FLAGS_traffic_matrix.empty());

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology), &node_order);
  builder.RemoveMultipleLinks();
  builder.ScaleCapacity(FLAGS_link_capacity_scale);
  nc::net::GraphStorage graph(builder);

  ctr::PcapTraceStore trace_store(FLAGS_pcap_trace_store);
  std::vector<ctr::BinSequence> all_bin_sequences;
  for (ctr::PcapDataTrace* trace : trace_store.AllTraces()) {
    all_bin_sequences.emplace_back(trace->ToSequence(trace->AllSlices()));
  }
  ctr::BinSequenceGenerator bin_sequence_generator(all_bin_sequences, 1000);

  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
      nc::lp::DemandMatrix::LoadRepetitaOrDie(
          nc::File::ReadFileToStringOrDie(FLAGS_traffic_matrix), node_order,
          &graph);
  demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);

  std::map<ctr::AggregateId, ctr::BinSequence> initial_sequences;
  for (const auto& matrix_element : demand_matrix->elements()) {
    ctr::BinSequence bin_sequence =
        ctr::BinsAtRate(matrix_element.demand,
                        std::chrono::milliseconds(FLAGS_period_duration_ms),
                        &bin_sequence_generator);
    LOG(ERROR) << bin_sequence.bin_count();

    initial_sequences.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(matrix_element.src, matrix_element.dst),
        std::forward_as_tuple(bin_sequence));
    break;
  }

  ctr::PathProvider path_provider(&graph);
  std::unique_ptr<ctr::Optimizer> opt;
  if (FLAGS_opt == "CTR") {
    opt = nc::make_unique<ctr::CTROptimizer>(&path_provider, false);
  } else if (FLAGS_opt == "B4") {
    opt = nc::make_unique<ctr::B4Optimizer>(&path_provider, false);
  } else if (FLAGS_opt == "B4(P)") {
    opt = nc::make_unique<ctr::B4Optimizer>(&path_provider, true);
  } else if (FLAGS_opt == "MinMax") {
    opt = nc::make_unique<ctr::MinMaxOptimizer>(&path_provider);
  } else if (FLAGS_opt == "SP") {
    opt = nc::make_unique<ctr::ShortestPathOptimizer>(&path_provider);
  }

  ctr::MeanScaleEstimatorFactory estimator_factory(
      {1.1, FLAGS_decay_factor, FLAGS_decay_factor, 10});
  ctr::RoutingSystem routing_system({}, opt.get(), &estimator_factory);

  nc::net::IPAddress controller_ip(100);
  nc::net::IPAddress device_ip_base(2000);
  nc::net::IPAddress tldr_ip_base(3000);
  nc::net::IPAddress flow_group_ip_base(4000);
  nc::net::DevicePortNumber enter_port(4000);
  nc::net::DevicePortNumber exit_port(5000);

  std::chrono::milliseconds round_duration(FLAGS_period_duration_ms);
  std::chrono::milliseconds poll_period(FLAGS_history_bin_size_ms);

  ctr::controller::Controller controller(controller_ip, &routing_system,
                                         &event_queue, &graph);
  ctr::MockSimDeviceFactory device_factory(enter_port, &event_queue);

  ctr::controller::NetworkContainerConfig containter_config(
      device_ip_base, tldr_ip_base, flow_group_ip_base, enter_port, exit_port,
      std::chrono::milliseconds(100), false);
  ctr::TLDRConfig tldr_config({}, nc::htsim::kWildIPAddress,
                              nc::htsim::kWildIPAddress, controller_ip,
                              round_duration, poll_period, 100);
  ctr::controller::NetworkContainer network_container(
      containter_config, tldr_config, &graph, &controller, &event_queue);

  // Will first add aggregates and populate the device factory.
  for (const auto& id_and_bin_sequence : initial_sequences) {
    const ctr::AggregateId& id = id_and_bin_sequence.first;
    const ctr::BinSequence& bin_sequence = id_and_bin_sequence.second;
    ctr::BinSequence from_start = bin_sequence.CutFromStart(round_duration);
    ctr::AggregateHistory init_history =
        from_start.GenerateHistory(poll_period);

    nc::htsim::MatchRuleKey key_for_aggregate =
        network_container.AddAggregate(id, init_history);

    // Will add the bin sequence to the source device.
    const std::string& src_device_id = graph.GetNode(id.src())->id();
    device_factory.AddBinSequence(src_device_id, key_for_aggregate,
                                  bin_sequence);

    // Will also add the reverse aggregate for ACKs.
    auto ack_history =
        GetDummyHistory(nc::net::Bandwidth::FromKBitsPerSecond(100), 10);
    network_container.AddAggregate(id.Reverse(), ack_history);
  }

  // Now that the device factory is ready, we can initialize the container.
  network_container.AddElementsFromGraph(&device_factory);

  // With devices and queues in place we can add flow groups to the container.
  for (const auto& id_and_bin_sequence : initial_sequences) {
    ctr::controller::TCPFlowGroup flow_group(
        10, std::chrono::milliseconds(20), std::chrono::milliseconds(50), 10000,
        std::chrono::milliseconds(100));
    flow_group.set_mean_object_size_fixed(true);
    flow_group.set_mean_wait_time_fixed(false);
    flow_group.AddKeyFrame(std::chrono::milliseconds(0), 10,
                           nc::net::Bandwidth::FromGBitsPerSecond(10));
    network_container.AddTCPFlowGroup(id_and_bin_sequence.first, flow_group);
  }
  network_container.InitAggregatesInController();

  ctr::NetInstrument net_instrument(
      network_container.queues(), network_container.flow_group_tcp_sources(),
      std::chrono::milliseconds(10), &event_queue);
  device_factory.Init();

  ProgressIndicator progress_indicator(std::chrono::milliseconds(100),
                                       &event_queue);
  event_queue.RunAndStopIn(std::chrono::seconds(90));
  nc::metrics::DefaultMetricManager()->PersistAllMetrics();
}
