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
#include "opt_eval.h"
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

static auto* path_stretch_ms =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint32_t, std::string, std::string, std::string>(
            "opt_path_stretch_ms",
            "How far away from the shortest path a path is (absolute)",
            "Topology", "Traffic matrix", "Optimizer");

static auto* path_sp_delay_ms =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint32_t, std::string, std::string, std::string>(
            "opt_path_sp_delay_ms", "Delay of the shortest path", "Topology",
            "Traffic matrix", "Optimizer");

static auto* path_stretch_rel =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string>(
            "opt_path_stretch_rel",
            "How far away from the shortest path a path is (relative)",
            "Topology", "Traffic matrix", "Optimizer");

static auto* path_flow_count =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint32_t, std::string, std::string, std::string>(
            "opt_path_flow_count", "Number of flows on each path", "Topology",
            "Traffic matrix", "Optimizer");

static auto* aggregate_path_count =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint32_t, std::string, std::string, std::string>(
            "opt_path_count", "Number of paths in each aggregate", "Topology",
            "Traffic matrix", "Optimizer");

static auto* link_utilization =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string>(
            "opt_link_utilization", "Per-link utilization", "Topology",
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

namespace ctr {

static constexpr size_t kTopAggregateFlowCount = 10000;

static std::unique_ptr<TrafficMatrix> FromDemandMatrix(
    const nc::lp::DemandMatrix& demand_matrix) {
  nc::net::Bandwidth max_demand = nc::net::Bandwidth::Zero();
  for (const auto& element : demand_matrix.elements()) {
    max_demand = std::max(max_demand, element.demand);
  }

  // Will assign flow counts so that the max bandwidth aggregate has
  // kTopAggregateFlowCount, and all other aggregates proportionally less.
  std::map<AggregateId, DemandAndFlowCount> demands_and_counts;
  for (const auto& element : demand_matrix.elements()) {
    size_t flow_count = kTopAggregateFlowCount * (element.demand / max_demand);
    flow_count = std::max(1ul, flow_count);

    AggregateId id(element.src, element.dst);
    demands_and_counts[id] = {element.demand, flow_count};
  }

  return nc::make_unique<TrafficMatrix>(demand_matrix.graph(),
                                        demands_and_counts);
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

  auto* path_stretch_handle = path_stretch_ms->GetHandle(topology, tm, opt);
  auto* path_sp_delay_handle = path_sp_delay_ms->GetHandle(topology, tm, opt);
  auto* path_flow_count_handle = path_flow_count->GetHandle(topology, tm, opt);
  auto* path_stretch_rel_handle =
      path_stretch_rel->GetHandle(topology, tm, opt);
  auto* unmet_demand_handle = path_unmet_demand->GetHandle(topology, tm, opt);
  auto* path_count_handle = aggregate_path_count->GetHandle(topology, tm, opt);

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

    microseconds shortest_path_delay = aggregate_id.GetSPDelay(*graph);
    milliseconds sp_delay_ms = duration_cast<milliseconds>(shortest_path_delay);

    // Will limit the delay at 1ms.
    sp_delay_ms = std::max(sp_delay_ms, milliseconds(1));

    size_t path_count = 0;
    for (const auto& route : routes) {
      const nc::net::Walk* path = route.first;
      double fraction = route.second;
      CHECK(fraction > 0);

      microseconds path_delay = path->delay();
      CHECK(path_delay >= shortest_path_delay);
      ++path_count;

      milliseconds path_delay_ms = duration_cast<milliseconds>(path_delay);
      path_delay_ms = std::max(path_delay_ms, milliseconds(1));

      nc::net::Bandwidth per_flow_rate =
          nc::FindOrDieNoPrint(per_flow_rates, path);
      nc::net::Bandwidth unmet = std::max(nc::net::Bandwidth::Zero(),
                                          required_per_flow - per_flow_rate);
      double delta_rel =
          static_cast<double>(path_delay_ms.count()) / sp_delay_ms.count();
      milliseconds delta = path_delay_ms - sp_delay_ms;

      path_sp_delay_handle->AddValue(path_delay_ms.count());
      path_stretch_handle->AddValue(delta.count());
      path_flow_count_handle->AddValue(fraction * total_num_flows);
      path_stretch_rel_handle->AddValue(delta_rel);
      unmet_demand_handle->AddValue(unmet.bps());

      for (nc::net::GraphLinkIndex link : path->links()) {
        link_to_total_load[link] += total_aggregate_demand * fraction;
      }
    }

    path_count_handle->AddValue(path_count);
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

static void RunOptimizers(const OptEvalInput& input) {
  const std::string& top_file = input.topology_file;
  const std::string& tm_file = input.tm_file;
  const nc::net::GraphStorage* graph = input.demand_matrix->graph();
  LOG(ERROR) << "Running " << top_file << " " << tm_file;

  // A separate PathProvider instance for B4. Want to run CTR with a fresh
  // PathProvider, but have to run B4 first.
  PathProvider b4_path_provider(graph);

  std::unique_ptr<TrafficMatrix> tm = FromDemandMatrix(*input.demand_matrix);
  std::unique_ptr<RoutingConfiguration> routing;
  routing = RunB4(*tm, top_file, tm_file, &b4_path_provider);

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

  std::vector<std::unique_ptr<nc::net::GraphStorage>> graphs;
  std::vector<ctr::OptEvalInput> to_process;
  std::tie(graphs, to_process) = ctr::GetOptEvalInputs();

  nc::RunInParallel<ctr::OptEvalInput>(
      to_process, [](const ctr::OptEvalInput& input) {
        ctr::RunOptimizers(input);
      }, FLAGS_threads);
}