#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <string>
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
#include "ctr.h"
#include "opt.h"
#include "opt_compare.h"
#include "path_provider.h"

DEFINE_string(topology_files, "", "A list of topology files");
DEFINE_string(matrix_files, "", "A list of matrix files");
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

static std::unique_ptr<ctr::TrafficMatrix> ScaleFlowCountsProportionaltely(
    const ctr::TrafficMatrix& from, const ctr::TrafficMatrix& to) {
  std::map<ctr::AggregateId, ctr::DemandAndFlowCount> new_demands;
  for (const auto& aggregate_and_demands : from.demands()) {
    const ctr::AggregateId& aggregate = aggregate_and_demands.first;
    const ctr::DemandAndFlowCount& from_demands = aggregate_and_demands.second;
    const ctr::DemandAndFlowCount& to_demands =
        nc::FindOrDieNoPrint(to.demands(), aggregate);

    // Will update the flow count in to_demands to preserve the
    // demand/flow_count ration from from_demands.
    size_t new_count = to_demands.first.Mbps() * from_demands.second /
                       from_demands.first.Mbps();
    new_demands[aggregate] = {to_demands.first, new_count};
  }

  return nc::make_unique<ctr::TrafficMatrix>(to.graph(), new_demands);
}

// Runs B4. If B4 is unable to fit the traffic will scale the matrix down to the
// point where it can.
static std::pair<std::unique_ptr<ctr::RoutingConfiguration>,
                 std::unique_ptr<ctr::RoutingConfiguration>>
RunB4(const ctr::TrafficMatrix& tm, ctr::PathProvider* path_provider,
      std::mt19937* rnd, double* scale_factor) {
  ctr::B4Optimizer b4_optimizer(path_provider, false);
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

    *scale_factor = scale;
    return {std::move(before), std::move(after)};
  }

  LOG(FATAL) << "Should not happen";
  return {};
}

static void PrintDelta(const ctr::RoutingConfiguration& before,
                       const ctr::RoutingConfiguration& after) {
  ctr::RoutingConfigurationDelta delta = before.GetDifference(after);
  for (const auto& aggregate_and_delta : delta.aggregates) {
    const ctr::AggregateId& aggregate = aggregate_and_delta.first;
    const ctr::AggregateDelta& aggregate_delta = aggregate_and_delta.second;

    const ctr::DemandAndFlowCount& demands_before =
        nc::FindOrDieNoPrint(before.demands(), aggregate);
    const ctr::DemandAndFlowCount& demands_after =
        nc::FindOrDieNoPrint(after.demands(), aggregate);

    if (aggregate_delta.fraction_delta > 0) {
      LOG(ERROR) << "Delta " << aggregate_delta.fraction_delta << " demand "
                 << demands_before.first.Mbps() << " -> "
                 << demands_after.first.Mbps() << " flow count "
                 << demands_before.second << " -> " << demands_after.second;
      LOG(ERROR) << before.AggregateToString(aggregate);
      LOG(ERROR) << after.AggregateToString(aggregate);
    }
  }
}

static void ParseMatrix(const ctr::TrafficMatrix& tm, size_t seed,
                        ctr::PathProvider* path_provider,
                        ctr::RoutingConfigDeltaInfo* b4_delta_info,
                        ctr::RoutingConfigDeltaInfo* ctr_delta_info,
                        ctr::RoutingConfigDeltaInfo* ctr_delta_p_info,
                        ctr::RoutingConfigDeltaInfo* ctr_delta_h_info,
                        double* scale_factor, double* b4_delay_delta,
                        double* h_delay_delta) {
  std::mt19937 rnd(seed);
  auto b4_before_and_after = RunB4(tm, path_provider, &rnd, scale_factor);
  ctr::RoutingConfigurationDelta b4_delta =
      b4_before_and_after.first->GetDifference(*b4_before_and_after.second);

  const ctr::TrafficMatrix* tm_before = b4_before_and_after.first.get();
  const ctr::TrafficMatrix* tm_after = b4_before_and_after.second.get();
  std::unique_ptr<ctr::TrafficMatrix> tm_proportionately_scaled =
      ScaleFlowCountsProportionaltely(*tm_before, *tm_after);

  // Need to also generate before/after for CTR using the same TMs as B4.
  ctr::CTROptimizer ctr_optimizer(path_provider);
  auto ctr_before = ctr_optimizer.Optimize(*tm_before);

  ctr::CTROptimizer ctr_optimizer_two(path_provider);
  auto ctr_after = ctr_optimizer_two.Optimize(*tm_after);
  ctr::RoutingConfigurationDelta ctr_delta =
      ctr_before->GetDifference(*ctr_after);
  PrintDelta(*ctr_before, *ctr_after);
  LOG(ERROR) << "\n\n";

  // Will run with the proportionately scaled matrix.
  ctr::CTROptimizer ctr_optimizer_three(path_provider);
  auto ctr_after_p = ctr_optimizer_three.Optimize(*tm_proportionately_scaled);
  ctr::RoutingConfigurationDelta ctr_delta_p =
      ctr_before->GetDifference(*ctr_after_p);

  //  nc::viz::HtmlPage page_before;
  //  ctr_before->ToHTML(&page_before);
  //  nc::File::WriteStringToFile(page_before.Construct(), "before.html");
  //
  //  nc::viz::HtmlPage page_after;
  //  ctr_after_p->ToHTML(&page_after);
  //  nc::File::WriteStringToFile(page_after.Construct(), "after.html");

  // Will also run with a heuristic.
  auto ctr_after_h = ctr_optimizer.OptimizeWithPrevious(*tm_after, *ctr_before);
  ctr::RoutingConfigurationDelta ctr_delta_h =
      ctr_before->GetDifference(*ctr_after_h);
  LOG(ERROR) << "Heuristic:";
  PrintDelta(*ctr_before, *ctr_after_h);

  nc::net::Delay best_delay = ctr_after->TotalPerFlowDelay();
  nc::net::Delay b4_delay = b4_before_and_after.second->TotalPerFlowDelay();
  nc::net::Delay h_delay = ctr_after_h->TotalPerFlowDelay();

  *b4_delay_delta =
      (b4_delay - best_delay).count() / static_cast<double>(best_delay.count());
  *h_delay_delta =
      (h_delay - best_delay).count() / static_cast<double>(best_delay.count());

  //  LOG(ERROR) << ctr_after_h->ToString();
  //  LOG(ERROR) << ctr_after->ToString();

  //  CHECK(b4_delay >= best_delay) << b4_delay.count() << " vs "
  //                                << best_delay.count();
  //
  CHECK(h_delay >= best_delay) << h_delay.count() << " vs "
                               << best_delay.count();

  b4_delta_info->Add(b4_delta);
  ctr_delta_info->Add(ctr_delta);
  ctr_delta_p_info->Add(ctr_delta_p);
  ctr_delta_h_info->Add(ctr_delta_h);
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
  if (!FLAGS_matrix_files.empty()) {
    return nc::Glob(FLAGS_matrix_files);
  }

  std::string matrix_location =
      nc::StringReplace(topology_file, ".graph", ".*.demands", true);
  return nc::Glob(matrix_location);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
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

    ctr::RoutingConfigDeltaInfo b4_delta_info;
    ctr::RoutingConfigDeltaInfo ctr_delta_info;
    ctr::RoutingConfigDeltaInfo ctr_delta_p_info;
    ctr::RoutingConfigDeltaInfo ctr_delta_h_info;

    // Extra info for each run.
    nc::viz::NpyArray extra_info(
        {{"topology", nc::viz::NpyArray::STRING},
         {"tm", nc::viz::NpyArray::STRING},
         {"seed", nc::viz::NpyArray::UINT64},
         {"downscale_factor", nc::viz::NpyArray::DOUBLE},
         {"b4_total_delay_delta", nc::viz::NpyArray::DOUBLE},
         {"ctr_total_delay_delta", nc::viz::NpyArray::DOUBLE}});

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
        double scale_factor = 1.0;
        double b4_total_delay_delta = 0.0;
        double ctr_h_total_delay_delta = 0.0;

        ParseMatrix(*scaled_tm, seed, &path_provider, &b4_delta_info,
                    &ctr_delta_info, &ctr_delta_p_info, &ctr_delta_h_info,
                    &scale_factor, &b4_total_delay_delta,
                    &ctr_h_total_delay_delta);
        extra_info.AddRow({topology_file, matrix_file, seed, scale_factor,
                           b4_total_delay_delta, ctr_h_total_delay_delta});
      }
    }

    b4_delta_info.AddPrefixToFieldNames("b4_");
    ctr_delta_info.AddPrefixToFieldNames("ctr_");
    ctr_delta_h_info.AddPrefixToFieldNames("ctrh_");
    ctr_delta_p_info.AddPrefixToFieldNames("ctrp_");

    extra_info.Combine(b4_delta_info)
        .Combine(ctr_delta_info)
        .Combine(ctr_delta_h_info)
        .Combine(ctr_delta_p_info);
    extra_info.ToDisk("out", topology_file != topology_files.front());
  }
}
