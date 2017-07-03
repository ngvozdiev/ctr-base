#include <gflags/gflags.h>
#include <chrono>
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
#include "common.h"
#include "controller.h"
#include "mean_est/mean_est.h"
#include "metrics/metrics.h"
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
DEFINE_uint64(initial_window_ms, 60000,
              "Initial window to guess a trace's rate");
DEFINE_double(decay_factor, 0.0, "How quickly to decay prediction");
DEFINE_double(link_capacity_scale, 1.0, "By how much to scale all links");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_string(opt, "CTR", "The optimizer to use");

static nc::htsim::MatchRuleKey KeyFromAggregateId(const ctr::AggregateId& id) {
  using namespace nc::htsim;
  using namespace nc::net;
  FiveTuple five_tuple(IPAddress(id.src()), IPAddress(id.dst()), kWildIPProto,
                       kWildAccessLayerPort, kWildAccessLayerPort);
  return {kWildPacketTag, kWildDevicePortNumber, {five_tuple}};
}

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
                        std::chrono::milliseconds(FLAGS_initial_window_ms),
                        &bin_sequence_generator);
    LOG(ERROR) << bin_sequence.bin_count();

    initial_sequences.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(matrix_element.src, matrix_element.dst),
        std::forward_as_tuple(bin_sequence));
    //    break;
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
  }

  ctr::MeanScaleEstimatorFactory estimator_factory(
      {1.1, FLAGS_decay_factor, FLAGS_decay_factor, 10});
  ctr::RoutingSystem routing_system({}, opt.get(), &estimator_factory);

  nc::net::IPAddress controller_ip(100);
  nc::net::IPAddress device_ip_base(2000);
  nc::net::IPAddress tldr_ip_base(3000);
  nc::net::DevicePortNumber enter_port(4000);
  nc::net::DevicePortNumber exit_port(5000);

  std::chrono::milliseconds round_duration(60000);
  std::chrono::milliseconds poll_period(100);

  ctr::controller::Controller controller(controller_ip, &routing_system,
                                         &event_queue, &graph);
  ctr::MockSimDeviceFactory device_factory(enter_port, &event_queue);
  for (const auto& id_and_bin_sequence : initial_sequences) {
    const ctr::AggregateId& id = id_and_bin_sequence.first;
    const ctr::BinSequence& bin_sequence = id_and_bin_sequence.second;
    nc::htsim::MatchRuleKey key_for_aggregate = KeyFromAggregateId(id);

    // Will add the bin sequence to the source device.
    const std::string& src_device_id = graph.GetNode(id.src())->id();
    device_factory.AddBinSequence(src_device_id, key_for_aggregate,
                                  bin_sequence);
  }

  ctr::controller::NetworkContainerConfig containter_config(
      device_ip_base, tldr_ip_base, enter_port, exit_port,
      std::chrono::milliseconds(100), false);
  ctr::TLDRConfig tldr_config({}, nc::htsim::kWildIPAddress,
                              nc::htsim::kWildIPAddress, controller_ip,
                              round_duration, poll_period, 100);
  ctr::controller::NetworkContainer network_container(
      containter_config, tldr_config, &graph, &controller, &device_factory,
      &event_queue);

  for (const auto& id_and_bin_sequence : initial_sequences) {
    const ctr::AggregateId& id = id_and_bin_sequence.first;
    const ctr::BinSequence& bin_sequence = id_and_bin_sequence.second;

    ctr::BinSequence from_start = bin_sequence.CutFromStart(round_duration);
    ctr::AggregateHistory init_history =
        from_start.GenerateHistory(poll_period);
    nc::htsim::MatchRuleKey key_for_aggregate = KeyFromAggregateId(id);

    ctr::controller::AggregateState aggregate_state(
        id.src(), id.dst(), key_for_aggregate, init_history);
    network_container.AddAggregate(id, aggregate_state);
  }

  ctr::NetInstrument net_instrument(network_container.queues(),
                                    std::chrono::milliseconds(1), &event_queue);
  network_container.Init();
  device_factory.Init();
  event_queue.RunAndStopIn(std::chrono::seconds(900));
}
