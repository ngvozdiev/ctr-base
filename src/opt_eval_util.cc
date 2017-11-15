#include <gflags/gflags.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/thread_runner.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/perfect_hash.h"
#include "demand_matrix_input.h"
#include "common.h"
#include "metrics/metrics.h"
#include "opt/ctr.h"
#include "opt/opt.h"
#include "opt/oversubscription_model.h"
#include "opt/path_provider.h"

using namespace std::chrono;

DEFINE_bool(scale_to_make_b4_fit, true,
            "If true will always scale down the TM to make B4 fit before "
            "running the rest of the optimizers.");
DEFINE_uint64(threads, 4, "Number of parallel threads to run");

static auto* aggregate_sp_delay_ms =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint32_t, std::string, std::string>(
            "aggregate_sp_delay_ms", "Delay of each aggregate's shortest path",
            "Topology", "Traffic matrix");

static auto* aggregate_rate_Mbps =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string>(
            "aggregate_rate_Mbps", "Rate of each aggregate", "Topology",
            "Traffic matrix");

static auto* aggregate_path_count =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint32_t, std::string, std::string, std::string>(
            "opt_path_count", "Number of paths in each aggregate", "Topology",
            "Traffic matrix", "Optimizer");

static auto* path_delay_ms =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint32_t, std::string, std::string, std::string>(
            "opt_path_delay_ms", "Delay of each path", "Topology",
            "Traffic matrix", "Optimizer");

static auto* path_flow_count =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint32_t, std::string, std::string, std::string>(
            "opt_path_flow_count", "Number of flows on each path", "Topology",
            "Traffic matrix", "Optimizer");

static auto* path_unmet_demand =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint64_t, std::string, std::string, std::string>(
            "opt_path_unmet_demand_bps",
            "How many bps of a single flow's demand are not satisfied",
            "Topology", "Traffic matrix", "Optimizer");

static auto* ctr_runtime_ms =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint64_t, std::string, std::string>(
            "ctr_runtime_ms", "How long it took CTR to run", "Topology",
            "Traffic matrix");

static auto* ctr_runtime_cached_ms =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint64_t, std::string, std::string>(
            "ctr_runtime_cached_ms", "How long it took CTR to run (cached)",
            "Topology", "Traffic matrix");

static auto* tm_scale_factor =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string>(
            "tm_scale_factor",
            "By how much the TM had to be scaled to make B4 fit", "Topology",
            "Traffic matrix");

static auto* total_delay_fraction =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string>(
            "total_delay_fraction", "Fraction of total delay at no headroom",
            "Topology", "Traffic matrix");

static auto* link_utilization =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string>(
            "opt_link_utilization", "Per-link utilization", "Topology",
            "Traffic matrix", "Optimizer");

namespace ctr {

static constexpr double kHeadroomStrideSize = 0.02;

static void RecordHeadroomVsDelay(const std::string& top_file,
                                  const std::string& tm_file,
                                  const TrafficMatrix& tm) {
  using namespace ctr;
  using namespace std::chrono;
  const nc::net::GraphStorage* graph = tm.graph();

  std::vector<double> values;
  for (double link_multiplier = 1.0; link_multiplier > 0.0;
       link_multiplier -= kHeadroomStrideSize) {
    if (!tm.ToDemandMatrix()->IsFeasible({}, link_multiplier)) {
      break;
    }

    PathProvider path_provider(graph);
    CTROptimizer ctr_optimizer(&path_provider, link_multiplier, false, false);
    std::unique_ptr<RoutingConfiguration> routing = ctr_optimizer.Optimize(tm);

    double max_utilization = routing->MaxLinkUtilization();
    if (max_utilization > link_multiplier + 0.001) {
      break;
    }

    double value = duration_cast<seconds>(routing->TotalPerFlowDelay()).count();
    values.emplace_back(value);
  }

  // The lowest delay is achieved at link_multiplier=1.0. Will normalize
  // everything by that.
  CHECK(!values.empty());
  double lowest_delay = values.front();
  auto* metric_handle = total_delay_fraction->GetHandle(top_file, tm_file);
  for (auto value : values) {
    metric_handle->AddValue(value / lowest_delay);
  }
}

static void RecordTrafficMatrixStats(const std::string& topology,
                                     const std::string& tm,
                                     const TrafficMatrix& traffic_matrix) {
  const nc::net::GraphStorage* graph = traffic_matrix.graph();
  auto* sp_delay_handle = aggregate_sp_delay_ms->GetHandle(topology, tm);
  auto* rate_handle = aggregate_rate_Mbps->GetHandle(topology, tm);
  for (const auto& aggregate_and_demand : traffic_matrix.demands()) {
    const AggregateId& aggregate_id = aggregate_and_demand.first;
    const DemandAndFlowCount& demand_and_flow_count =
        aggregate_and_demand.second;

    microseconds shortest_path_delay = aggregate_id.GetSPDelay(*graph);
    milliseconds sp_delay_ms = duration_cast<milliseconds>(shortest_path_delay);

    // Will limit the delay at 1ms.
    sp_delay_ms = std::max(sp_delay_ms, milliseconds(1));
    sp_delay_handle->AddValue(sp_delay_ms.count());
    rate_handle->AddValue(demand_and_flow_count.first.Mbps());
  }
}

static void RecordRoutingConfig(const std::string& topology,
                                const std::string& tm, const std::string& opt,
                                const RoutingConfiguration& routing) {
  using namespace std::chrono;
  const nc::net::GraphStorage* graph = routing.graph();

  // A map from a link to the total load over the link.
  std::map<nc::net::GraphLinkIndex, nc::net::Bandwidth> link_to_total_load;

  OverSubModel model(routing);
  const std::map<const nc::net::Walk*, nc::net::Bandwidth>& per_flow_rates =
      model.per_flow_bandwidth_map();

  auto* path_flow_count_handle = path_flow_count->GetHandle(topology, tm, opt);
  auto* unmet_demand_handle = path_unmet_demand->GetHandle(topology, tm, opt);
  auto* path_count_handle = aggregate_path_count->GetHandle(topology, tm, opt);
  auto* path_delay = path_delay_ms->GetHandle(topology, tm, opt);

  for (const auto& aggregate_and_aggregate_output : routing.routes()) {
    const AggregateId& aggregate_id = aggregate_and_aggregate_output.first;
    const std::vector<RouteAndFraction>& routes =
        aggregate_and_aggregate_output.second;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(routing.demands(), aggregate_id);

    size_t total_num_flows = demand_and_flow_count.second;
    nc::net::Bandwidth total_aggregate_demand = demand_and_flow_count.first;
    nc::net::Bandwidth required_per_flow =
        total_aggregate_demand / total_num_flows;

    for (const auto& route : routes) {
      const nc::net::Walk* path = route.first;
      double fraction = route.second;
      CHECK(fraction > 0);

      microseconds delay = path->delay();
      milliseconds delay_ms = duration_cast<milliseconds>(delay);
      delay_ms = std::max(delay_ms, milliseconds(1));
      path_delay->AddValue(delay_ms.count());

      nc::net::Bandwidth per_flow_rate =
          nc::FindOrDieNoPrint(per_flow_rates, path);
      nc::net::Bandwidth unmet = std::max(nc::net::Bandwidth::Zero(),
                                          required_per_flow - per_flow_rate);

      path_flow_count_handle->AddValue(fraction * total_num_flows);
      unmet_demand_handle->AddValue(unmet.bps());

      for (nc::net::GraphLinkIndex link : path->links()) {
        link_to_total_load[link] += total_aggregate_demand * fraction;
      }
    }

    path_count_handle->AddValue(routes.size());
  }

  auto* link_load_handle = link_utilization->GetHandle(topology, tm, opt);
  nc::net::GraphLinkSet links_seen;
  for (const auto& link_and_total_load : link_to_total_load) {
    nc::net::GraphLinkIndex link_index = link_and_total_load.first;
    links_seen.Insert(link_index);
    const nc::net::GraphLink* link = graph->GetLink(link_index);

    nc::net::Bandwidth total_load = link_and_total_load.second;
    link_load_handle->AddValue(total_load / link->bandwidth());
  }

  for (nc::net::GraphLinkIndex link_index : graph->AllLinks()) {
    if (!links_seen.Contains(link_index)) {
      link_load_handle->AddValue(0);
    }
  }
}

static bool Fits(const ctr::RoutingConfiguration& routing) {
  const std::set<ctr::AggregateId>& aggregates_no_fit =
      ctr::OverSubModel(routing).aggregates_no_fit();
  return aggregates_no_fit.empty();
}

// Runs B4. If B4 is unable to fit the traffic will scale the matrix down to the
// point where it can.
static std::unique_ptr<ctr::RoutingConfiguration> RunB4(
    const ctr::TrafficMatrix& tm, const std::string& topology_string,
    const std::string& tm_string, ctr::PathProvider* path_provider) {
  ctr::B4Optimizer b4_optimizer(path_provider, false, 1.0);
  double scale = 1.01;
  while (scale > 0) {
    // Will scale it down by 1%.
    scale -= 0.01;

    auto scaled_tm = tm.ScaleDemands(scale, {});
    // B4 should run fine on both the scaled and the randomized TMs.
    auto routing = b4_optimizer.Optimize(*scaled_tm);
    if (FLAGS_scale_to_make_b4_fit && !Fits(*routing)) {
      continue;
    }

    tm_scale_factor->GetHandle(topology_string, tm_string)->AddValue(scale);
    return routing;
  }

  LOG(FATAL) << "Should not happen";
  return {};
}

static void RunOptimizers(const DemandMatrixAndFilename& input) {
  const std::string& top_file = input.topology_file;
  const std::string& tm_file = input.file;
  const nc::net::GraphStorage* graph = input.demand_matrix->graph();
  LOG(ERROR) << "Running " << top_file << " " << tm_file;

  // A separate PathProvider instance for B4. Want to run CTR with a fresh
  // PathProvider, but have to run B4 first.
  PathProvider b4_path_provider(graph);

  std::unique_ptr<TrafficMatrix> tm =
      TrafficMatrix::ProportionalFromDemandMatrix(*input.demand_matrix);
  std::unique_ptr<RoutingConfiguration> routing;
  routing = RunB4(*tm, top_file, tm_file, &b4_path_provider);

  RecordTrafficMatrixStats(top_file, tm_file, *tm);
  RecordHeadroomVsDelay(top_file, tm_file, *tm);
  RecordRoutingConfig(top_file, tm_file, "B4", *routing);

  PathProvider path_provider(graph);
  CTROptimizer ctr_optimizer(&path_provider, 1.0, false, false);
  CTROptimizer ctr_optimizer_no_flow_counts(&path_provider, 1.0, false, true);
  MinMaxOptimizer minmax_optimizer(&path_provider, 1.0, false);
  MinMaxOptimizer minmax_low_delay_optimizer(&path_provider, 1.0, true);
  MinMaxPathBasedOptimizer minmax_ksp_optimizer(&path_provider, 1.0, true, 10);
  B4Optimizer b4_flow_count_optimizer(&path_provider, true, 1.0);

  auto ctr_start = high_resolution_clock::now();
  routing = ctr_optimizer.Optimize(*tm);
  auto ctr_duration = high_resolution_clock::now() - ctr_start;

  RecordRoutingConfig(top_file, tm_file, "CTR", *routing);

  ctr_start = high_resolution_clock::now();
  ctr_optimizer.Optimize(*tm);
  auto ctr_cached_duration = high_resolution_clock::now() - ctr_start;

  routing = ctr_optimizer_no_flow_counts.Optimize(*tm);
  RecordRoutingConfig(top_file, tm_file, "CTRNFC", *routing);

  routing = minmax_optimizer.Optimize(*tm);
  RecordRoutingConfig(top_file, tm_file, "MinMax", *routing);

  routing = minmax_low_delay_optimizer.Optimize(*tm);
  RecordRoutingConfig(top_file, tm_file, "MinMaxLD", *routing);

  routing = minmax_ksp_optimizer.Optimize(*tm);
  RecordRoutingConfig(top_file, tm_file, "MinMaxK10", *routing);

  routing = b4_flow_count_optimizer.Optimize(*tm);
  RecordRoutingConfig(top_file, tm_file, "B4FC", *routing);

  auto* handle = ctr_runtime_ms->GetHandle(top_file, tm_file);
  handle->AddValue(duration_cast<milliseconds>(ctr_duration).count());
  handle = ctr_runtime_cached_ms->GetHandle(top_file, tm_file);
  handle->AddValue(duration_cast<milliseconds>(ctr_cached_duration).count());
}

}  // namespace ctr

using TopologyAndMatrix = std::tuple<std::string, std::string, double>;
int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::metrics::InitMetrics();

  // Will switch off timestamps.
  auto timestamp_provider =
      ::nc::make_unique<nc::metrics::NullTimestampProvider>();
  nc::metrics::DefaultMetricManager()->set_timestamp_provider(
      std::move(timestamp_provider));

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> to_process;
  std::tie(topologies, to_process) = ctr::GetDemandMatrixInputs();

  nc::RunInParallel<ctr::DemandMatrixAndFilename>(
      to_process, [](const ctr::DemandMatrixAndFilename& input) {
        ctr::RunOptimizers(input);
      }, FLAGS_threads);
}
