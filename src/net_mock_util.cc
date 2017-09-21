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

using namespace std::chrono;

DEFINE_string(topology, "", "A file with a topology");
DEFINE_string(traffic_matrix, "", "A file with a traffic matrix");
DEFINE_string(pcap_trace_store, "trace_store.pb",
              "A file with information about .pcap traces");
DEFINE_string(pcap_trace_fit_store, "",
              "A file with a series of PBTraceToFitRate protobufs");
DEFINE_uint64(period_duration_ms, 60000, "Length of the period");
DEFINE_uint64(history_bin_size_ms, 100, "How big each history bin is");
DEFINE_double(decay_factor, 0.0, "How quickly to decay prediction");
DEFINE_double(estimator_headroom, 1.1, "Fixed headroom for the estimator");
DEFINE_double(link_capacity_scale, 1.0,
              "By how much to scale all links' bandwidth");
DEFINE_double(link_delay_scale, 1.3, "By how much to scale all links' delay");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_uint64(duration_ms, 90000, "For how long to run (in simulated time)");
DEFINE_double(exceed_probability, 0.001,
              "Probability convolution to exceed rate");

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

static uint64_t FlowCountFromBinsSequence(const ctr::BinSequence& bin_sequence,
                                          ctr::PcapDataBinCache* bin_cache) {
  std::vector<uint64_t> syns;

  for (const auto& bin :
       bin_sequence.AccumulateBins(bin_sequence.bin_size(), bin_cache)) {
    syns.emplace_back(bin.flows_enter);
  }

  return ctr::GetFlowCountFromSyns(syns);
}

static void HandleDefault(
    const nc::net::GraphStorage& graph,
    std::map<ctr::AggregateId, std::unique_ptr<ctr::BinSequence>>&&
        initial_sequences,
    nc::net::DevicePortNumber enter_port, milliseconds poll_period,
    milliseconds round_duration, nc::EventQueue* event_queue,
    ctr::controller::NetworkContainer* network_container) {
  // Will first add aggregates and populate the device factory.
  ctr::MockSimDeviceFactory device_factory(enter_port, event_queue);
  for (auto& id_and_bin_sequence : initial_sequences) {
    const ctr::AggregateId& id = id_and_bin_sequence.first;
    std::unique_ptr<ctr::BinSequence>& bin_sequence =
        id_and_bin_sequence.second;
    std::unique_ptr<ctr::BinSequence> from_start =
        bin_sequence->CutFromStart(round_duration);

    uint64_t flow_count =
        FlowCountFromBinsSequence(*from_start, device_factory.bin_cache());
    ctr::AggregateHistory init_history = from_start->GenerateHistory(
        poll_period, flow_count, device_factory.bin_cache());

    nc::htsim::MatchRuleKey key_for_aggregate =
        network_container->AddAggregate(id, init_history);

    // Will add the bin sequence to the source device.
    const std::string& src_device_id = graph.GetNode(id.src())->id();
    device_factory.AddBinSequence(src_device_id, key_for_aggregate,
                                  std::move(bin_sequence));

    // Will also add the reverse aggregate for ACKs.
    //    auto ack_history =
    //        ctr::GetDummyHistory(nc::net::Bandwidth::FromKBitsPerSecond(100),
    //                             milliseconds(FLAGS_history_bin_size_ms),
    //                             milliseconds(FLAGS_period_duration_ms), 10);
    //    network_container->AddAggregate(id.Reverse(), ack_history);
  }

  // Records per-path stats.
  ctr::InputPacketObserver input_packet_observer(
      network_container->controller(), milliseconds(10), event_queue);

  // Monitors queues.
  ctr::OutputPacketObserver output_packet_observer(event_queue);

  // Now that the device factory is ready, we can initialize the container.
  network_container->AddElementsFromGraph(
      &device_factory, &input_packet_observer, &output_packet_observer);
  network_container->InitAggregatesInController();

  ctr::NetInstrument net_instrument(network_container->internal_queues(),
                                    network_container->flow_group_tcp_sources(),
                                    milliseconds(10), event_queue);
  device_factory.Init();
  event_queue->RunAndStopIn(milliseconds(FLAGS_duration_ms));
  output_packet_observer.RecordDist();
}

static std::string StripExtension(const std::string& filename) {
  std::size_t found = filename.find_last_of(".");
  return filename.substr(0, found);
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

  nc::net::Delay diameter = graph.Stats().sp_delay_percentiles.back();
  if (diameter < milliseconds(10)) {
    LOG(FATAL) << "Graph diameter too small ("
               << duration_cast<milliseconds>(diameter).count() << "ms)";
  }

  ctr::PcapTraceStore trace_store(FLAGS_pcap_trace_store);
  std::string fit_store = FLAGS_pcap_trace_fit_store;
  if (fit_store.empty()) {
    fit_store = nc::StrCat(StripExtension(FLAGS_traffic_matrix), "_fit.pb");
  }
  CHECK(nc::File::Exists(fit_store)) << " fit store does not exist "
                                     << fit_store;

  ctr::PcapTraceFitStore trace_fit_store(fit_store, &trace_store);
  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
      nc::lp::DemandMatrix::LoadRepetitaOrDie(
          nc::File::ReadFileToStringOrDie(FLAGS_traffic_matrix), node_order,
          &graph);
  demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);

  size_t i = -1;
  std::map<ctr::AggregateId, std::unique_ptr<ctr::BinSequence>>
      initial_sequences;
  for (const auto& matrix_element : demand_matrix->elements()) {
    LOG(INFO) << "Looking for traces to match " << matrix_element.demand.Mbps()
              << "Mbps";
    std::unique_ptr<ctr::BinSequence> bin_sequence =
        trace_fit_store.GetBinSequence(matrix_element.demand);
    CHECK(bin_sequence);

    std::unique_ptr<ctr::BinSequence> bin_sequence_extended =
        trace_store.ExtendBinSequence(*bin_sequence);

    initial_sequences.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(matrix_element.src, matrix_element.dst),
        std::forward_as_tuple(std::move(bin_sequence_extended)));
    LOG(INFO) << "Bin sequence " << ++i << " / "
              << demand_matrix->elements().size();
  }

  ctr::PathProvider path_provider(&graph);
  ctr::CTROptimizer opt(&path_provider, 0.95, true, false);

  ctr::ProbModelConfig prob_model_config;
  prob_model_config.exceed_probability = FLAGS_exceed_probability;

  ctr::MeanScaleEstimatorFactory estimator_factory(
      {FLAGS_estimator_headroom, FLAGS_decay_factor, FLAGS_decay_factor, 10});
  ctr::RoutingSystemConfig routing_system_config;
  routing_system_config.prob_model_config = prob_model_config;
  ctr::RoutingSystem routing_system(routing_system_config, &opt,
                                    &estimator_factory);

  nc::net::IPAddress controller_ip(100);
  nc::net::IPAddress device_ip_base(20000);
  nc::net::IPAddress tldr_ip_base(30000);
  nc::net::IPAddress flow_group_ip_base(40000);
  nc::net::DevicePortNumber enter_port(4000);
  nc::net::DevicePortNumber exit_port(5000);

  milliseconds round_duration(FLAGS_period_duration_ms);
  milliseconds poll_period(FLAGS_history_bin_size_ms);

  ctr::controller::Controller controller(controller_ip, &routing_system,
                                         &event_queue, &graph);
  ctr::controller::NetworkContainerConfig containter_config(
      device_ip_base, tldr_ip_base, flow_group_ip_base, enter_port, exit_port,
      milliseconds(100), milliseconds(1000000), false);

  nc::ThresholdEnforcerPolicy te_policy;
  ctr::TLDRConfig tldr_config(te_policy, prob_model_config,
                              nc::htsim::kWildIPAddress,
                              nc::htsim::kWildIPAddress, controller_ip,
                              round_duration, poll_period, 100);
  ctr::controller::NetworkContainer network_container(
      containter_config, tldr_config, &graph, &controller, &event_queue);

  nc::htsim::ProgressIndicator progress_indicator(milliseconds(100),
                                                  &event_queue);
  HandleDefault(graph, std::move(initial_sequences), enter_port, poll_period,
                round_duration, &event_queue, &network_container);

  nc::metrics::DefaultMetricManager()->PersistAllMetrics();
}
