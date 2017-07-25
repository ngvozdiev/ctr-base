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
#include "common.h"
#include "metrics/metrics.h"
#include "opt/ctr.h"
#include "opt/opt.h"
#include "opt/oversubscription_model.h"
#include "opt/path_provider.h"

DEFINE_string(topology_files, "", "Topology files");
DEFINE_double(link_capacity_scale, 1.0, "By how much to scale all links");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_double(delay_scale, 1.0, "By how much to scale the delays of all links");
DEFINE_string(output, "opt_eval_out", "Output directory");
DEFINE_uint64(threads, 4, "Number of parallel threads to run");

static auto* path_stretch_ms =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint32_t, std::string, std::string, std::string>(
            "opt_path_stretch_ms",
            "How far away from the shortest path a path is (absolute)",
            "Topology", "Traffic matrix", "Optimizer");

static auto* path_stretch_rel =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string, std::string>(
            "opt_path_stretch_rel",
            "How far away from the shortest path a path is (relative)",
            "Topology", "Traffic matrix", "Optimizer");

static auto* path_flow_count =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint32_t, std::string, std::string, std::string>(
            "opt_path_flow_count", "Number of flows on each path", "Topology",
            "Traffic matrix", "Optimizer");

static auto* aggregate_path_count =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint32_t, std::string, std::string, std::string>(
            "opt_path_count", "Number of paths in each aggregate", "Topology",
            "Traffic matrix", "Optimizer");

static auto* link_utilization =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string, std::string>(
            "opt_link_utilization", "Per-link utilization", "Topology",
            "Traffic matrix", "Optimizer");

static auto* path_unmet_demand =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string, std::string, std::string>(
            "opt_path_unmet_demand_bps",
            "How many bps of a single flow's demand are not satisfied",
            "Topology", "Traffic matrix", "Optimizer");

static auto* ctr_runtime_ms =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string, std::string>(
            "ctr_runtime_ms", "How long it took CTR to run", "Topology",
            "Traffic matrix");

static auto* ctr_runtime_cached_ms =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string, std::string>(
            "ctr_runtime_cached_ms", "How long it took CTR to run (cached)",
            "Topology", "Traffic matrix");

static auto* tm_scale_factor =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "tm_scale_factor",
            "By how much the TM had to be scaled to make B4 fit", "Topology",
            "Traffic matrix");

namespace ctr {

static std::unique_ptr<TrafficMatrix> FromDemandMatrix(
    const nc::lp::DemandMatrix& demand_matrix) {
  std::map<AggregateId, DemandAndFlowCount> demands_and_counts;
  for (const auto& element : demand_matrix.elements()) {
    AggregateId id(element.src, element.dst);
    demands_and_counts[id] = {element.demand, 1000ul};
  }

  return nc::make_unique<TrafficMatrix>(demand_matrix.graph(),
                                        demands_and_counts);
}

static std::mutex global_mutex;

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

      path_stretch_handle->AddValue(path_delay_ms.count());
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
    if (!Fits(*routing)) {
      continue;
    }

    {
      std::unique_lock<std::mutex> lock(global_mutex);
      tm_scale_factor->GetHandle(topology_string, tm_string)->AddValue(scale);
    }
    return routing;
  }

  LOG(FATAL) << "Should not happen";
  return {};
}

static void RunOptimizers(const std::string& topology_file,
                          const std::string& tm_file) {
  using namespace std::chrono;
  LOG(ERROR) << "Running " << topology_file << " " << tm_file;

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(topology_file), &node_order);
  builder.RemoveMultipleLinks();
  builder.ScaleCapacity(FLAGS_link_capacity_scale);
  builder.ScaleDelay(FLAGS_delay_scale);
  nc::net::GraphStorage graph(builder);

  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
      nc::lp::DemandMatrix::LoadRepetitaOrDie(
          nc::File::ReadFileToStringOrDie(tm_file), node_order, &graph);
  demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);

  std::string top_file_trimmed = nc::Split(topology_file, "/").back();
  std::string tm_file_trimmed = nc::Split(tm_file, "/").back();

  // A separate PathProvider instance for B4. Want to run CTR with a fresh
  // PathProvider, but have to run B4 first.
  PathProvider b4_path_provider(&graph);

  std::unique_ptr<TrafficMatrix> tm = FromDemandMatrix(*demand_matrix);
  std::unique_ptr<RoutingConfiguration> routing;
  routing = RunB4(*tm, top_file_trimmed, tm_file_trimmed, &b4_path_provider);

  {
    std::unique_lock<std::mutex> lock(global_mutex);
    RecordRoutingConfig(top_file_trimmed, tm_file_trimmed, "B4", *routing);
  }

  PathProvider path_provider(&graph);
  CTROptimizer ctr_optimizer(&path_provider, 1.0, false);
  MinMaxOptimizer minmax_optimizer(&path_provider, 1.0, false);
  MinMaxOptimizer minmax_low_delay_optimizer(&path_provider, 1.0, true);
  B4Optimizer b4_flow_count_optimizer(&path_provider, true, 1.0);

  auto ctr_start = high_resolution_clock::now();
  routing = ctr_optimizer.Optimize(*tm);
  auto ctr_duration = high_resolution_clock::now() - ctr_start;

  {
    std::unique_lock<std::mutex> lock(global_mutex);
    RecordRoutingConfig(top_file_trimmed, tm_file_trimmed, "CTR", *routing);
  }

  ctr_start = std::chrono::high_resolution_clock::now();
  ctr_optimizer.Optimize(*tm);
  auto ctr_cached_duration =
      std::chrono::high_resolution_clock::now() - ctr_start;

  routing = minmax_optimizer.Optimize(*tm);

  {
    std::unique_lock<std::mutex> lock(global_mutex);
    RecordRoutingConfig(top_file_trimmed, tm_file_trimmed, "MinMax", *routing);
  }

  routing = minmax_low_delay_optimizer.Optimize(*tm);

  {
    std::unique_lock<std::mutex> lock(global_mutex);
    RecordRoutingConfig(top_file_trimmed, tm_file_trimmed, "MinMaxLD",
                        *routing);
  }

  routing = b4_flow_count_optimizer.Optimize(*tm);

  std::unique_lock<std::mutex> lock(global_mutex);
  RecordRoutingConfig(top_file_trimmed, tm_file_trimmed, "B4FC", *routing);

  auto* handle = ctr_runtime_ms->GetHandle(top_file_trimmed, tm_file_trimmed);
  handle->AddValue(duration_cast<milliseconds>(ctr_duration).count());
  handle = ctr_runtime_cached_ms->GetHandle(top_file_trimmed, tm_file_trimmed);
  handle->AddValue(duration_cast<milliseconds>(ctr_cached_duration).count());
}

}  // namespace ctr

static std::vector<std::string> GetTopologyFiles() {
  std::vector<std::string> out;
  std::vector<std::string> split = nc::Split(FLAGS_topology_files, ",");
  for (const std::string& piece : split) {
    std::vector<std::string> files = nc::Glob(piece);
    out.insert(out.end(), files.begin(), files.end());
  }

  return out;
}

static std::vector<std::string> GetMatrixFiles(
    const std::string& topology_file) {
  std::string matrix_location =
      nc::StringReplace(topology_file, ".graph", ".*.demands", true);
  return nc::Glob(matrix_location);
}

using TopologyAndMatrix = std::pair<std::string, std::string>;
int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::metrics::InitMetrics();

  std::vector<std::string> topology_files = GetTopologyFiles();
  CHECK(!topology_files.empty());

  std::vector<TopologyAndMatrix> to_process;
  for (const std::string& topology_file : topology_files) {
    std::vector<std::string> matrix_files = GetMatrixFiles(topology_file);
    if (matrix_files.empty()) {
      LOG(ERROR) << "No matrices for " << topology_file;
      continue;
    }

    for (const std::string& matrix_file : matrix_files) {
      to_process.emplace_back(topology_file, matrix_file);
    }
  }

  nc::RunInParallel<TopologyAndMatrix>(
      to_process, [](const TopologyAndMatrix& tm_and_matrix) {
        ctr::RunOptimizers(tm_and_matrix.first, tm_and_matrix.second);
      }, FLAGS_threads);
}
