#include <gflags/gflags.h>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "ncode/common.h"
#include "ncode/event_queue.h"
#include "ncode/file.h"
#include "ncode/htsim/match.h"
#include "ncode/htsim/network.h"
#include "ncode/logging.h"
#include "ncode/lp/demand_matrix.h"
#include "ncode/net/net_common.h"
#include "ncode/net/net_gen.h"
#include "ncode/strutil.h"
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
#include "prob_model/dist_model.h"
#include "routing_system.h"
#include "tldr.h"

using namespace std::chrono;

DEFINE_string(topology, "", "A file with a topology");
DEFINE_string(traffic_matrix, "", "A file with a traffic matrix");
DEFINE_string(pcap_trace_store, "trace_store.pb",
              "A file with information about .pcap traces");
DEFINE_string(pcap_trace_fit_store, "",
              "A file with a series of PBTraceToFitRate protobufs");
DEFINE_uint64(period_duration_ms, 5000, "Length of the period");
DEFINE_uint64(periods_in_history, 12, "Periods in history");
DEFINE_uint64(history_bin_size_ms, 100, "How big each history bin is");
DEFINE_double(decay_factor, 0.0, "How quickly to decay prediction");
DEFINE_double(link_capacity_multiplier, 0.9,
              "Link capacity multiplier during optimization");
DEFINE_double(link_capacity_scale, 1.0,
              "By how much to scale all links' bandwidth");
DEFINE_double(link_delay_scale, 1.3, "By how much to scale all links' delay");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_uint64(duration_ms, 90000, "For how long to run (in simulated time)");
DEFINE_double(exceed_probability, 0.001,
              "Probability convolution to exceed rate");
DEFINE_bool(sp_opt, false, "If true will run SP routing instead of CTR");
DEFINE_double(demand_scale, 1.0, "Will scale all packet traces by this much");
DEFINE_bool(limited_opt, true,
            "Whether or not to attempt to avoid churn at the expense of a "
            "little optimality");
DEFINE_uint64(cache_prepopulate, 0,
              "If not 0 will pre-populate the cache with this many bins from "
              "each trace/slices combination");

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

static void HandleApproximate(
    std::map<ctr::AggregateId, ctr::BinSequence>&& initial_sequences,
    ctr::RoutingSystem* routing_system) {
  milliseconds round_duration(FLAGS_period_duration_ms);
  milliseconds poll_period(FLAGS_history_bin_size_ms);
  milliseconds total_duration(FLAGS_duration_ms);

  ctr::PcapDataBinCache cache;
  if (FLAGS_cache_prepopulate != 0) {
    std::map<const ctr::PcapDataTrace*, ctr::TraceSliceSet> to_cache;
    for (const auto& aggregate_and_init_sequence : initial_sequences) {
      const ctr::BinSequence& sequence = aggregate_and_init_sequence.second;
      for (const ctr::BinSequence::TraceAndSlice& trace_and_slice :
           sequence.traces()) {
        to_cache[trace_and_slice.trace].insert(trace_and_slice.slice);
      }
    }

    for (const auto& trace_and_slices : to_cache) {
      cache.Populate(trace_and_slices.first, trace_and_slices.second, 0,
                     FLAGS_cache_prepopulate);
    }
  }

  ctr::NetMock net_mock(std::move(initial_sequences), round_duration,
                        poll_period, total_duration, FLAGS_periods_in_history,
                        routing_system);
  net_mock.Run(&cache);
}

static std::string StripExtension(const std::string& filename) {
  std::size_t found = filename.find_last_of(".");
  return filename.substr(0, found);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::metrics::InitMetrics();

  CHECK(!FLAGS_topology.empty());
  CHECK(!FLAGS_traffic_matrix.empty());

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology), &node_order);
  builder.RemoveMultipleLinks();
  builder.ScaleCapacity(FLAGS_link_capacity_scale);
  builder.ScaleDelay(FLAGS_link_delay_scale);
  nc::net::GraphStorage graph(builder, node_order);

  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
      nc::lp::DemandMatrix::LoadRepetitaFileOrDie(FLAGS_traffic_matrix,
                                                  node_order, &graph);
  demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);

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

  size_t i = -1;
  std::map<ctr::AggregateId, ctr::BinSequence> initial_sequences;
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
        std::forward_as_tuple(bin_sequence_extended->traces()));
    LOG(INFO) << "Bin sequence " << ++i << " / "
              << demand_matrix->elements().size();
  }

  ctr::PathProvider path_provider(&graph);
  std::unique_ptr<ctr::Optimizer> opt;
  if (FLAGS_sp_opt) {
    opt = nc::make_unique<ctr::ShortestPathOptimizer>(&path_provider);
  } else {
    opt = nc::make_unique<ctr::CTROptimizer>(&path_provider,
                                             FLAGS_link_capacity_multiplier,
                                             FLAGS_limited_opt, false);
  }

  ctr::ProbModelConfig prob_model_config;
  prob_model_config.exceed_probability = FLAGS_exceed_probability;

  ctr::MeanScaleEstimatorFactory estimator_factory(
      {1.0, FLAGS_decay_factor, FLAGS_decay_factor, 10});
  ctr::RoutingSystemConfig routing_system_config;
  routing_system_config.prob_model_config = prob_model_config;
  routing_system_config.store_to_metrics = false;
  routing_system_config.link_capacity_multiplier =
      FLAGS_link_capacity_multiplier +
      (1 - FLAGS_link_capacity_multiplier) / 3.0;
  ctr::RoutingSystem routing_system(routing_system_config, opt.get(),
                                    &estimator_factory);

  HandleApproximate(std::move(initial_sequences), &routing_system);
  nc::metrics::DefaultMetricManager()->PersistAllMetrics();
}
