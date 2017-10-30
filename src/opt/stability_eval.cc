#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/lp/mc_flow.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/strutil.h"
#include "../common.h"
#include "../metrics/metrics.h"
#include "ctr.h"
#include "opt.h"
#include "oversubscription_model.h"
#include "path_provider.h"

static auto* stability_volume_delta =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string, std::string>(
            "stability_volume_delta", "Fraction of volume that changed",
            "Topology", "Traffic matrix", "Optimizer");

static auto* stability_flow_count_delta =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string, std::string>(
            "stability_flow_count_delta", "Fraction of volume that changed",
            "Topology", "Traffic matrix", "Optimizer");

static auto* stability_volume_delta_longer_path =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string, std::string>(
            "stability_volume_delta_longer_path",
            "Fraction of volume that moved to longer path", "Topology",
            "Traffic matrix", "Optimizer");

static auto* stability_flow_count_delta_longer_path =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string, std::string>(
            "stability_flow_count_delta_longer_path",
            "Fraction of flows that moved to longer path", "Topology",
            "Traffic matrix", "Optimizer");

static auto* stability_add_count =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint32_t, std::string, std::string, std::string>(
            "stability_add_count", "Number of ADDs", "Topology",
            "Traffic matrix", "Optimizer");

static auto* stability_update_count =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint32_t, std::string, std::string, std::string>(
            "stability_update_count", "Number of UPDATEs", "Topology",
            "Traffic matrix", "Optimizer");

static auto* stability_remove_count =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint32_t, std::string, std::string, std::string>(
            "stability_remove_count", "Number of REMOVEs", "Topology",
            "Traffic matrix", "Optimizer");

static auto* stability_scale_factor =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "stability_scale_factor",
            "By how much the TM had to be scaled to make B4 fit", "Topology",
            "Traffic matrix");

DEFINE_string(topology_files, "", "A list of topology files");
DEFINE_double(demand_fraction, 0.05, "By how much to vary demand");
DEFINE_double(flow_count_fraction, 0.05, "By how much to vary flow counts.");
DEFINE_uint64(aggregate_count, 1, "How many aggregates to change");
DEFINE_uint64(try_count, 1, "How many tries to perform");
DEFINE_uint64(seed, 1ul, "Seed to use for the RNG");
DEFINE_bool(scale_to_fit_b4, true, "Scales down matrices until B4 fits.");
DEFINE_uint64(max_matrix_size, 500, "Ignore matrices with > aggregates.");
DEFINE_double(matrix_scale, 1.0, "Matrix demands will be scaled by this much");

static bool Fits(const ctr::RoutingConfiguration& routing) {
  const std::set<ctr::AggregateId>& aggregates_no_fit =
      ctr::OverSubModel(routing).aggregates_no_fit();
  return aggregates_no_fit.empty();
}

void RecordDeltas(const std::string& topology, const std::string& tm,
                  const std::string& opt,
                  const ctr::RoutingConfigurationDelta& delta) {
  using namespace std::chrono;

  auto* volume_delta_handle =
      stability_volume_delta->GetHandle(topology, tm, opt);
  auto* flow_count_delta_handle =
      stability_flow_count_delta->GetHandle(topology, tm, opt);
  auto* volume_delta_lp_handle =
      stability_volume_delta_longer_path->GetHandle(topology, tm, opt);
  auto* flow_count_delta_lp_handle =
      stability_flow_count_delta_longer_path->GetHandle(topology, tm, opt);
  auto* add_count_handle = stability_add_count->GetHandle(topology, tm, opt);
  auto* update_count_handle =
      stability_update_count->GetHandle(topology, tm, opt);
  auto* remove_count_handle =
      stability_remove_count->GetHandle(topology, tm, opt);

  volume_delta_handle->AddValue(delta.total_volume_fraction_delta);
  volume_delta_lp_handle->AddValue(delta.total_volume_fraction_on_longer_path);
  flow_count_delta_handle->AddValue(delta.total_flow_fraction_delta);
  flow_count_delta_lp_handle->AddValue(
      delta.total_flow_fraction_on_longer_path);

  size_t route_adds;
  size_t route_removals;
  size_t route_updates;
  std::tie(route_adds, route_removals, route_updates) = delta.TotalRoutes();
  add_count_handle->AddValue(route_adds);
  update_count_handle->AddValue(route_updates);
  remove_count_handle->AddValue(route_removals);
}

// Runs B4. If B4 is unable to fit the traffic will scale the matrix down to the
// point where it can.
static std::pair<std::unique_ptr<ctr::RoutingConfiguration>,
                 std::unique_ptr<ctr::RoutingConfiguration>>
RunB4(const ctr::TrafficMatrix& tm, const std::string& topology_string,
      const std::string& tm_string, ctr::PathProvider* path_provider,
      std::mt19937* rnd) {
  ctr::B4Optimizer b4_optimizer(path_provider, false, 1.0);
  double scale = 1.01;
  while (scale > 0) {
    // Will scale it down by 1%.
    scale -= 0.01;

    auto scaled_tm = tm.ScaleDemands(scale, {});
    // B4 should run fine on both the scaled and the randomized TMs.
    auto before = b4_optimizer.Optimize(*scaled_tm);
    if (FLAGS_scale_to_fit_b4 && !Fits(*before)) {
      continue;
    }

    auto rnd_tm =
        scaled_tm->Randomize(FLAGS_demand_fraction, FLAGS_flow_count_fraction,
                             FLAGS_aggregate_count, rnd);
    auto after = b4_optimizer.Optimize(*rnd_tm);
    if (FLAGS_scale_to_fit_b4 && !Fits(*after)) {
      continue;
    }

    stability_scale_factor->GetHandle(topology_string, tm_string)
        ->AddValue(scale);
    return {std::move(before), std::move(after)};
  }

  LOG(FATAL) << "Should not happen";
  return {};
}

static void ParseMatrix(const ctr::TrafficMatrix& tm,
                        const std::string& topology_string,
                        const std::string& tm_string, size_t seed,
                        ctr::PathProvider* path_provider) {
  std::mt19937 rnd(seed);
  auto b4_before_and_after =
      RunB4(tm, topology_string, tm_string, path_provider, &rnd);
  ctr::RoutingConfigurationDelta b4_delta =
      b4_before_and_after.first->GetDifference(*b4_before_and_after.second);

  const ctr::TrafficMatrix* tm_before = b4_before_and_after.first.get();
  const ctr::TrafficMatrix* tm_after = b4_before_and_after.second.get();

  // Need to also generate before/after for CTR using the same TMs as B4.
  ctr::CTROptimizer ctr_optimizer(path_provider, 1.0, true, false);
  auto ctr_before = ctr_optimizer.Optimize(*tm_before);

  // Will also run with a heuristic.
  std::unique_ptr<ctr::RoutingConfiguration> ctr_after_with_limits;
  std::unique_ptr<ctr::RoutingConfiguration> ctr_after;
  std::tie(ctr_after_with_limits, ctr_after) =
      ctr_optimizer.OptimizeAndReturnUnlimitedRun(*tm_after);

  ctr::RoutingConfigurationDelta ctr_delta_limits =
      ctr_before->GetDifference(*ctr_after_with_limits);
  ctr::RoutingConfigurationDelta ctr_delta =
      ctr_before->GetDifference(*ctr_after);

  ctr::MinMaxOptimizer minmax_opt(path_provider, 1.0, true);
  auto minmax_before = minmax_opt.Optimize(*tm_before);
  auto minmax_after = minmax_opt.Optimize(*tm_after);
  ctr::RoutingConfigurationDelta minmax_delta =
      minmax_before->GetDifference(*minmax_after);

  ctr::MinMaxPathBasedOptimizer minmax_pb_opt(path_provider, 1.0, true, 10);
  auto minmax_pb_before = minmax_pb_opt.Optimize(*tm_before);
  auto minmax_pb_after = minmax_pb_opt.Optimize(*tm_after);
  ctr::RoutingConfigurationDelta minmax_pb_delta =
      minmax_pb_before->GetDifference(*minmax_pb_after);

  RecordDeltas(topology_string, tm_string, "MinMax(LD)", minmax_delta);
  RecordDeltas(topology_string, tm_string, "MinMax(K10)", minmax_pb_delta);
  RecordDeltas(topology_string, tm_string, "CTR", ctr_delta);
  RecordDeltas(topology_string, tm_string, "CTR_LIM", ctr_delta_limits);
  RecordDeltas(topology_string, tm_string, "B4", b4_delta);
}

// The TM that we load will have no flow counts. Need some out of thin air.
static std::map<nc::lp::SrcAndDst, size_t> GetFlowCountMap(
    const nc::lp::DemandMatrix& demand_matrix) {
  std::map<nc::lp::SrcAndDst, size_t> out;
  for (const auto& element : demand_matrix.elements()) {
    // Will assign flow counts in a way that yields the same per-flow bandwidth.
    double count = element.demand.Mbps();
    out[{element.src, element.dst}] = std::max(1.0, count);
  }

  return out;
}

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

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::metrics::InitMetrics();

  std::vector<std::string> topology_files = GetTopologyFiles();
  CHECK(!topology_files.empty());

  for (const std::string& topology_file : topology_files) {
    std::vector<std::string> nodes_in_order;
    nc::net::GraphBuilder graph_builder = nc::net::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(topology_file), &nodes_in_order);
    graph_builder.RemoveMultipleLinks();
    nc::net::GraphStorage graph(graph_builder);

    std::vector<std::string> matrix_files = GetMatrixFiles(topology_file);
    if (matrix_files.empty()) {
      LOG(ERROR) << "No matrices for " << topology_file;
      continue;
    }

    for (const std::string& matrix_file : matrix_files) {
      LOG(INFO) << "Processing " << topology_file << " : " << matrix_file;
      auto demand_matrix = nc::lp::DemandMatrix::LoadRepetitaOrDie(
          nc::File::ReadFileToStringOrDie(matrix_file), nodes_in_order, &graph);
      ctr::TrafficMatrix traffic_matrix(*demand_matrix,
                                        GetFlowCountMap(*demand_matrix));
      if (traffic_matrix.demands().size() > FLAGS_max_matrix_size) {
        continue;
      }

      // Will scale the TM here.
      auto scaled_tm = traffic_matrix.ScaleDemands(FLAGS_matrix_scale, {});

      ctr::PathProvider path_provider(&graph);
      for (size_t i = 0; i < FLAGS_try_count; ++i) {
        size_t seed = FLAGS_seed + i;
        ParseMatrix(*scaled_tm, nc::File::ExtractFileName(matrix_file),
                    nc::File::ExtractFileName(topology_file), seed,
                    &path_provider);
      }
    }
  }
}
