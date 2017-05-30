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
DEFINE_double(demand_fraction, 0.05, "By how much to vary demand");
DEFINE_double(flow_count_fraction, std::numeric_limits<double>::max(),
              "By how much to vary flow counts. If set to max will update flow "
              "counts to be proportional to demand");
DEFINE_uint64(aggregate_count, 1, "How many aggregates to change");
DEFINE_uint64(try_count, 1, "How many tries to perform");
DEFINE_uint64(seed, 1ul, "Seed to use for the RNG");

static bool Fits(const ctr::RoutingConfiguration& routing) {
  const std::set<ctr::AggregateId>& aggregates_no_fit =
      ctr::OverSubModel(routing).aggregates_no_fit();
  return aggregates_no_fit.empty();
}

// Runs B4. If B4 is unable to fit the traffic will scale the matrix down to the
// point where it can.
static std::pair<std::unique_ptr<ctr::RoutingConfiguration>,
                 std::unique_ptr<ctr::RoutingConfiguration>>
RunB4(const ctr::TrafficMatrix& tm, ctr::PathProvider* path_provider,
      std::mt19937* rnd, double* scale_factor) {
  ctr::B4Optimizer b4_optimizer(path_provider);
  //  std::unique_ptr<ctr::TrafficMatrix> scaled_tm =
  //      nc::make_unique<ctr::TrafficMatrix>(tm.graph(), tm.demands());
  //  while (true) {
  //    auto before = b4_optimizer.Optimize(*scaled_tm);
  //    std::set<ctr::AggregateId> aggregates_no_fit =
  //        ctr::OverSubModel(*before).aggregates_no_fit();
  //    if (!aggregates_no_fit.empty()) {
  //      scaled_tm = scaled_tm->ScaleDemands(0.95, aggregates_no_fit);
  //      continue;
  //    }
  //
  //    auto rnd_tm =
  //        scaled_tm->Randomize(FLAGS_demand_fraction,
  //        FLAGS_flow_count_fraction,
  //                             FLAGS_aggregate_count, rnd);
  //    auto after = b4_optimizer.Optimize(*rnd_tm);
  //    aggregates_no_fit = ctr::OverSubModel(*after).aggregates_no_fit();
  //    if (!aggregates_no_fit.empty()) {
  //      scaled_tm = scaled_tm->ScaleDemands(0.95, aggregates_no_fit);
  //      continue;
  //    }
  //
  //    return {std::move(before), std::move(after)};
  //  }

  double scale = 1.01;
  while (scale > 0) {
    // Will scale it down by 1%.
    scale -= 0.01;

    auto scaled_tm = tm.ScaleDemands(scale, {});
    // B4 should run fine on both the scaled and the randomized TMs.
    auto before = b4_optimizer.Optimize(*scaled_tm);
    if (!Fits(*before)) {
      continue;
    }

    auto rnd_tm =
        scaled_tm->Randomize(FLAGS_demand_fraction, FLAGS_flow_count_fraction,
                             FLAGS_aggregate_count, rnd);
    auto after = b4_optimizer.Optimize(*rnd_tm);
    if (!Fits(*after)) {
      continue;
    }

    *scale_factor = scale;
    return {std::move(before), std::move(after)};
  }

  LOG(FATAL) << "Should not happen";
  return {};
}

static void ParseMatrix(const ctr::TrafficMatrix& tm, size_t seed,
                        ctr::PathProvider* path_provider,
                        ctr::RoutingConfigDeltaInfo* b4_delta_info,
                        ctr::RoutingConfigDeltaInfo* ctr_delta_info,
                        ctr::RoutingConfigDeltaInfo* ctr_delta_h_info,
                        double* scale_factor) {
  std::mt19937 rnd(seed);
  auto b4_before_and_after = RunB4(tm, path_provider, &rnd, scale_factor);
  ctr::RoutingConfigurationDelta b4_delta =
      b4_before_and_after.first->GetDifference(*b4_before_and_after.second);

  // Need to also generate before/after for CTR using the same TMs as B4.
  ctr::CTROptimizer ctr_optimizer(path_provider);
  auto ctr_before = ctr_optimizer.Optimize(*b4_before_and_after.first);
  auto ctr_after = ctr_optimizer.Optimize(*b4_before_and_after.second);
  ctr::RoutingConfigurationDelta ctr_delta =
      ctr_before->GetDifference(*ctr_after);

  nc::viz::HtmlPage page1;
  ctr_before->ToHTML(&page1);
  nc::File::WriteStringToFile(page1.Construct(), "before.html");

  nc::viz::HtmlPage page2;
  ctr_after->ToHTML(&page2);
  nc::File::WriteStringToFile(page2.Construct(), "after.html");

  // Will also run with a heuristic.
  auto ctr_after_h = ctr_optimizer.OptimizeWithPrevious(
      *b4_before_and_after.second, *ctr_before);
  ctr::RoutingConfigurationDelta ctr_delta_h =
      ctr_before->GetDifference(*ctr_after_h);

  b4_delta_info->Add(b4_delta);
  ctr_delta_info->Add(ctr_delta);
  ctr_delta_h_info->Add(ctr_delta_h);
}

// The TM that we load will have no flow counts. Need some out of thin air.
static std::map<nc::lp::SrcAndDst, size_t> GetFlowCountMap(
    const nc::lp::DemandMatrix& demand_matrix) {
  std::map<nc::lp::SrcAndDst, size_t> out;
  for (const auto& element : demand_matrix.elements()) {
    out[{element.src, element.dst}] = 1000;
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
  std::vector<std::string> topology_files = GetTopologyFiles();
  CHECK(!topology_files.empty());

  ctr::RoutingConfigDeltaInfo b4_delta_info;
  ctr::RoutingConfigDeltaInfo ctr_delta_info;
  ctr::RoutingConfigDeltaInfo ctr_delta_h_info;

  // Extra info for each run.
  nc::viz::NpyArray extra_info(
      {{"topology", nc::viz::NpyArray::STRING},
       {"tm", nc::viz::NpyArray::STRING},
       {"seed", nc::viz::NpyArray::UINT64},
       {"downscale_factor", nc::viz::NpyArray::DOUBLE}});

  for (const std::string& topology_file : topology_files) {
    std::vector<std::string> nodes_in_order;
    nc::net::GraphBuilder graph_builder = nc::net::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(topology_file), &nodes_in_order);
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

      ctr::PathProvider path_provider(&graph);
      for (size_t i = 0; i < FLAGS_try_count; ++i) {
        size_t seed = FLAGS_seed + i;
        double scale_factor;
        ParseMatrix(traffic_matrix, seed, &path_provider, &b4_delta_info,
                    &ctr_delta_info, &ctr_delta_h_info, &scale_factor);
        extra_info.AddRow({topology_file, matrix_file, seed, scale_factor});
      }
    }
  }

  auto combined = nc::viz::NpyArray::Combine(
      nc::viz::NpyArray::Combine(
          nc::viz::NpyArray::Combine(
              extra_info, b4_delta_info.AddPrefixToFieldNames("b4_")),
          ctr_delta_info.AddPrefixToFieldNames("ctr_")),
      ctr_delta_h_info.AddPrefixToFieldNames("ctrh_"));
  combined.ToDisk("out");
}
