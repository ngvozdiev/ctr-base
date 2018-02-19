#include <gflags/gflags.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/net/net_common.h>
#include <stddef.h>
#include <stdint.h>
#include <iostream>
#include <memory>
#include <random>
#include <tuple>
#include <vector>

#include "../common.h"
#include "../demand_matrix_input.h"
#include "../metrics/metrics.h"
#include "../topology_input.h"
#include "ctr.h"
#include "opt.h"
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

static auto* stability_total_flow_delay_delta =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string, std::string>(
            "stability_total_flow_delay_delta",
            "Change in the total sum of delays experienced by all flows.",
            "Topology", "Traffic matrix", "Optimizer");

using namespace ctr;

DEFINE_double(demand_fraction, 0.05, "By how much to vary demand");
DEFINE_double(flow_count_fraction, 0.05, "By how much to vary flow counts.");
DEFINE_uint64(aggregate_count, 1, "How many aggregates to change");
DEFINE_uint64(try_count, 1, "How many tries to perform");
DEFINE_uint64(stability_eval_seed, 1ul, "Seed to use for the RNG");

DEFINE_bool(skip_trivial, true, "Skips trivially satisfiable TMs");

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
  auto* total_delay_delta_handle =
      stability_total_flow_delay_delta->GetHandle(topology, tm, opt);

  volume_delta_handle->AddValue(delta.total_volume_fraction_delta);
  volume_delta_lp_handle->AddValue(delta.total_volume_fraction_on_longer_path);
  flow_count_delta_handle->AddValue(delta.total_flow_fraction_delta);
  flow_count_delta_lp_handle->AddValue(
      delta.total_flow_fraction_on_longer_path);
  total_delay_delta_handle->AddValue(delta.total_per_flow_delay_delta);

  size_t route_adds;
  size_t route_removals;
  size_t route_updates;
  std::tie(route_adds, route_removals, route_updates) = delta.TotalRoutes();
  add_count_handle->AddValue(route_adds);
  update_count_handle->AddValue(route_updates);
  remove_count_handle->AddValue(route_removals);
}

static void ParseMatrix(
    const ctr::DemandMatrixAndFilename* demand_matrix_and_filename,
    size_t seed) {
  std::mt19937 rnd(seed);
  const std::string& top_file = demand_matrix_and_filename->topology_file;
  const std::string& tm_file = demand_matrix_and_filename->file;
  nc::lp::DemandMatrix* demand_matrix =
      demand_matrix_and_filename->demand_matrix.get();

  const nc::net::GraphStorage* graph = demand_matrix->graph();
  LOG(ERROR) << "Running " << top_file << " " << tm_file;

  std::unique_ptr<ctr::TrafficMatrix> tm =
      ctr::TrafficMatrix::DistributeFromDemandMatrix(*demand_matrix);
  auto tm_after =
      tm->Randomize(FLAGS_demand_fraction, FLAGS_flow_count_fraction,
                    FLAGS_aggregate_count, &rnd);

  PathProvider path_provider(graph);
  CTROptimizer ctr_optimizer_limited(&path_provider, 1.0, true, false);
  MinMaxOptimizer minmax_low_delay_optimizer(&path_provider, 1.0, true);
  MinMaxPathBasedOptimizer minmax_ksp_optimizer(&path_provider, 1.0, true, 10);
  B4Optimizer b4_optimizer(&path_provider, false, 1.0);

  auto ctr_before = ctr_optimizer_limited.Optimize(*tm);
  std::unique_ptr<ctr::RoutingConfiguration> ctr_after_with_limits;
  std::unique_ptr<ctr::RoutingConfiguration> ctr_after;
  std::tie(ctr_after_with_limits, ctr_after) =
      ctr_optimizer_limited.OptimizeAndReturnUnlimitedRun(*tm_after);

  ctr::RoutingConfigurationDelta ctr_delta_limits =
      ctr_before->GetDifference(*ctr_after_with_limits);
  ctr::RoutingConfigurationDelta ctr_delta =
      ctr_before->GetDifference(*ctr_after);

  auto minmax_before = minmax_low_delay_optimizer.Optimize(*tm);
  auto minmax_after = minmax_low_delay_optimizer.Optimize(*tm_after);
  ctr::RoutingConfigurationDelta minmax_delta =
      minmax_before->GetDifference(*minmax_after);

  auto minmax_pb_before = minmax_ksp_optimizer.Optimize(*tm);
  auto minmax_pb_after = minmax_ksp_optimizer.Optimize(*tm_after);
  ctr::RoutingConfigurationDelta minmax_pb_delta =
      minmax_pb_before->GetDifference(*minmax_pb_after);

  auto b4_before = b4_optimizer.Optimize(*tm);
  auto b4_after = b4_optimizer.Optimize(*tm_after);
  ctr::RoutingConfigurationDelta b4_delta = b4_before->GetDifference(*b4_after);

  RecordDeltas(top_file, tm_file, "MinMaxLD", minmax_delta);
  RecordDeltas(top_file, tm_file, "MinMaxK10", minmax_pb_delta);
  RecordDeltas(top_file, tm_file, "CTR", ctr_delta);
  RecordDeltas(top_file, tm_file, "CTR_LIM", ctr_delta_limits);
  RecordDeltas(top_file, tm_file, "B4", b4_delta);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::metrics::InitMetrics();

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::tie(topologies, matrices) =
      ctr::GetDemandMatrixInputs(FLAGS_skip_trivial);

  for (const ctr::DemandMatrixAndFilename& input : matrices) {
    for (size_t i = 0; i < FLAGS_try_count; ++i) {
      size_t seed = FLAGS_stability_eval_seed + i;
      ParseMatrix(&input, seed);
    }
  }
}
