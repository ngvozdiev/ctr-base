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
DEFINE_double(link_capacity_scale, 1.0,
              "By how much to scale all links' bandwidth");
DEFINE_double(link_delay_scale, 1.3, "By how much to scale all links' delay");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_string(opt, "CTR", "The optimizer to use");
DEFINE_bool(quick, false, "Whether to perform packet-level simulation or not");
DEFINE_bool(disable_rto, false, "If true will disable TCP's rto mechanism");
DEFINE_uint64(duration_sec, 90, "For how long to run (in simulated time)");
DEFINE_uint64(
    tcp_initial_cwnd_pkts, 4,
    "How many packets should there be in the inital congestion window");
DEFINE_uint64(object_size_bytes, 10000, "How many byte each object should be");
DEFINE_bool(simulate_initial_handshake, true,
            "Whether or not to simulate an initial handshake");
DEFINE_uint64(tcp_tracer_flow_max_count, 50,
              "How many tracer TCP flows to add to the largest-volume "
              "aggregate. Other aggregates' flow count is proportional.");
DEFINE_bool(
    disable_fast_optimization_requests, false,
    "If true will disable fast optimization requests to the controller.");
DEFINE_bool(pin_mean, false,
            "If true will disable the controller's probability model-based "
            "scaling of aggregates before optimization. All aggregates will "
            "take their mean level.");
DEFINE_bool(pin_max, false,
            "If true will disable the controller's probability model-based "
            "scaling of aggregates before optimization. All aggregates will "
            "take their max level.");

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

static std::map<ctr::AggregateId, size_t> GetTCPTracerFlowCounts(
    const nc::lp::DemandMatrix& demand_matrix) {
  nc::net::Bandwidth max_demand = nc::net::Bandwidth::Zero();
  for (const auto& element : demand_matrix.elements()) {
    max_demand = std::max(max_demand, element.demand);
  }
  CHECK(max_demand != nc::net::Bandwidth::Zero());

  std::map<ctr::AggregateId, size_t> out;
  for (const auto& element : demand_matrix.elements()) {
    double f = element.demand.Mbps() / max_demand.Mbps();
    size_t tcp_flow_count = f * FLAGS_tcp_tracer_flow_max_count;
    tcp_flow_count = std::max(1ul, tcp_flow_count);
    out[{element.src, element.dst}] = tcp_flow_count;
  }

  return out;
}

static uint64_t FlowCountFromBinsSequence(ctr::BinSequence* bin_sequence) {
  std::vector<uint64_t> syns;

  for (const auto& bin :
       bin_sequence->AccumulateBins(bin_sequence->bin_size())) {
    syns.emplace_back(bin.flows_enter);
  }

  return ctr::GetFlowCountFromSyns(syns);
}

static void HandleDefault(
    const std::map<ctr::AggregateId, ctr::BinSequence>& initial_sequences,
    const std::map<ctr::AggregateId, size_t>& tracer_flow_counts,
    const nc::net::GraphStorage& graph, nc::net::DevicePortNumber enter_port,
    std::chrono::milliseconds poll_period,
    std::chrono::milliseconds round_duration, nc::EventQueue* event_queue,
    ctr::controller::NetworkContainer* network_container) {
  // Will first add aggregates and populate the device factory.
  ctr::MockSimDeviceFactory device_factory(enter_port, event_queue);
  for (const auto& id_and_bin_sequence : initial_sequences) {
    const ctr::AggregateId& id = id_and_bin_sequence.first;
    const ctr::BinSequence& bin_sequence = id_and_bin_sequence.second;
    ctr::BinSequence from_start = bin_sequence.CutFromStart(round_duration);

    uint64_t flow_count = FlowCountFromBinsSequence(&from_start);
    ctr::AggregateHistory init_history =
        from_start.GenerateHistory(poll_period, flow_count);

    nc::htsim::MatchRuleKey key_for_aggregate =
        network_container->AddAggregate(id, init_history);

    // Will add the bin sequence to the source device.
    const std::string& src_device_id = graph.GetNode(id.src())->id();
    device_factory.AddBinSequence(src_device_id, key_for_aggregate,
                                  bin_sequence);

    // Will also add the reverse aggregate for ACKs.
    auto ack_history = ctr::GetDummyHistory(
        nc::net::Bandwidth::FromKBitsPerSecond(100),
        std::chrono::milliseconds(FLAGS_history_bin_size_ms),
        std::chrono::milliseconds(FLAGS_period_duration_ms), 10);
    network_container->AddAggregate(id.Reverse(), ack_history);
  }

  // Records per-path stats.
  ctr::InputPacketObserver packet_observer(network_container->controller(),
                                           std::chrono::milliseconds(10),
                                           event_queue);

  // Now that the device factory is ready, we can initialize the container.
  network_container->AddElementsFromGraph(&device_factory, &packet_observer);

  // With devices and queues in place we can add flow groups to the container.
  for (const auto& id_and_bin_sequence : initial_sequences) {
    const ctr::AggregateId& id = id_and_bin_sequence.first;
    size_t count = nc::FindOrDieNoPrint(tracer_flow_counts, id);

    ctr::controller::TCPFlowGroup flow_group(
        count, std::chrono::milliseconds(10), std::chrono::milliseconds(20),
        FLAGS_object_size_bytes, std::chrono::milliseconds(200));
    flow_group.set_mean_object_size_fixed(true);
    flow_group.set_mean_wait_time_fixed(false);
    flow_group.AddKeyFrame(std::chrono::milliseconds(0), 10,
                           nc::net::Bandwidth::FromGBitsPerSecond(10));
    network_container->AddTCPFlowGroup(id, flow_group);
  }
  network_container->InitAggregatesInController();

  ctr::NetInstrument net_instrument(network_container->internal_queues(),
                                    network_container->flow_group_tcp_sources(),
                                    std::chrono::milliseconds(10), event_queue);
  device_factory.Init();
  event_queue->RunAndStopIn(std::chrono::seconds(FLAGS_duration_sec));
}

static void HandleQuick(
    const std::map<ctr::AggregateId, ctr::BinSequence>& initial_sequences,
    const nc::net::GraphStorage& graph, std::chrono::milliseconds poll_period,
    std::chrono::milliseconds round_duration, nc::EventQueue* event_queue,
    ctr::controller::NetworkContainer* network_container) {
  // Will first add aggregates and populate the device factory.
  ctr::MockDeviceFactory device_factory(network_container->controller());
  for (const auto& id_and_bin_sequence : initial_sequences) {
    const ctr::AggregateId& id = id_and_bin_sequence.first;
    const ctr::BinSequence& bin_sequence = id_and_bin_sequence.second;
    ctr::BinSequence from_start = bin_sequence.CutFromStart(round_duration);

    uint64_t flow_count = FlowCountFromBinsSequence(&from_start);
    ctr::AggregateHistory init_history =
        from_start.GenerateHistory(poll_period, flow_count);

    nc::htsim::MatchRuleKey key_for_aggregate =
        network_container->AddAggregate(id, init_history);

    // Will add the bin sequence to the source device.
    const std::string& src_device_id = graph.GetNode(id.src())->id();
    device_factory.AddBinSequence(src_device_id, key_for_aggregate,
                                  bin_sequence);
  }

  // Records per-path stats.
  ctr::InputPacketObserver packet_observer(network_container->controller(),
                                           std::chrono::milliseconds(10),
                                           event_queue);

  // Now that the device factory is ready, we can initialize the container.
  network_container->AddElementsFromGraph(&device_factory, &packet_observer);
  network_container->InitAggregatesInController();

  ctr::NetInstrument net_instrument(network_container->internal_queues(),
                                    network_container->flow_group_tcp_sources(),
                                    std::chrono::milliseconds(10), event_queue);
  event_queue->RunAndStopIn(std::chrono::seconds(FLAGS_duration_sec));
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
  builder.ScaleDelay(FLAGS_link_delay_scale);
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
    opt = nc::make_unique<ctr::CTROptimizer>(&path_provider, 0.95, true);
  } else if (FLAGS_opt == "B4") {
    opt = nc::make_unique<ctr::B4Optimizer>(&path_provider, false, 0.95);
  } else if (FLAGS_opt == "B4(P)") {
    opt = nc::make_unique<ctr::B4Optimizer>(&path_provider, true, 0.95);
  } else if (FLAGS_opt == "MinMax") {
    opt = nc::make_unique<ctr::MinMaxOptimizer>(&path_provider, 0.9);
  } else if (FLAGS_opt == "SP") {
    opt = nc::make_unique<ctr::ShortestPathOptimizer>(&path_provider);
  }

  ctr::MeanScaleEstimatorFactory estimator_factory(
      {1.05, FLAGS_decay_factor, FLAGS_decay_factor, 10});
  ctr::RoutingSystemConfig routing_system_config;
  routing_system_config.pin_max = FLAGS_pin_max;
  routing_system_config.pin_mean = FLAGS_pin_mean;
  ctr::RoutingSystem routing_system(routing_system_config, opt.get(),
                                    &estimator_factory);

  nc::net::IPAddress controller_ip(100);
  nc::net::IPAddress device_ip_base(20000);
  nc::net::IPAddress tldr_ip_base(30000);
  nc::net::IPAddress flow_group_ip_base(40000);
  nc::net::DevicePortNumber enter_port(4000);
  nc::net::DevicePortNumber exit_port(5000);

  std::chrono::milliseconds tcp_rto_timer_period(10);
  if (FLAGS_disable_rto) {
    tcp_rto_timer_period = std::chrono::hours(9999);
  }

  std::chrono::milliseconds round_duration(FLAGS_period_duration_ms);
  std::chrono::milliseconds poll_period(FLAGS_history_bin_size_ms);

  ctr::controller::Controller controller(controller_ip, &routing_system,
                                         &event_queue, &graph);
  ctr::controller::NetworkContainerConfig containter_config(
      device_ip_base, tldr_ip_base, flow_group_ip_base, enter_port, exit_port,
      std::chrono::milliseconds(100), tcp_rto_timer_period, false);
  containter_config.tcp_config.inital_cwnd_size = FLAGS_tcp_initial_cwnd_pkts;
  containter_config.tcp_config.simulate_initial_handshake =
      FLAGS_simulate_initial_handshake;

  nc::ThresholdEnforcerPolicy te_policy;
  ctr::TLDRConfig tldr_config(te_policy, {}, nc::htsim::kWildIPAddress,
                              nc::htsim::kWildIPAddress, controller_ip,
                              round_duration, poll_period, 100,
                              FLAGS_disable_fast_optimization_requests);
  ctr::controller::NetworkContainer network_container(
      containter_config, tldr_config, &graph, &controller, &event_queue);

  nc::htsim::ProgressIndicator progress_indicator(
      std::chrono::milliseconds(100), &event_queue);
  if (FLAGS_quick) {
    HandleQuick(initial_sequences, graph, poll_period, round_duration,
                &event_queue, &network_container);
  } else {
    HandleDefault(initial_sequences, GetTCPTracerFlowCounts(*demand_matrix),
                  graph, enter_port, poll_period, round_duration, &event_queue,
                  &network_container);
  }
  nc::metrics::DefaultMetricManager()->PersistAllMetrics();
}
