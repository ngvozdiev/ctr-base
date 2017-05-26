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

DEFINE_string(topology_file, "", "A topology file");
DEFINE_string(matrices, "", "One or more matrices");
DEFINE_double(demand_fraction, 0.05, "By how much to vary demand");
DEFINE_double(flow_count_fraction, std::numeric_limits<double>::max(),
              "By how much to vary flow counts. If set to max will update flow "
              "counts to be proportional to demand");
DEFINE_uint64(aggregate_count, 1, "How many aggregates to change");
DEFINE_uint64(try_count, 1, "How many tries to perform");
DEFINE_string(opt, "CTR", "Optimizer to run. One of CTR,MinMax,B4");
DEFINE_uint64(seed, 1ul, "Seed to use for the RNG");

template <typename T>
static std::string PercentilesToString(
    std::vector<T>* values,
    const std::vector<size_t>& percenile_indices = {50, 90, 95, 100}) {
  std::vector<T> percentiles = nc::Percentiles(values);
  std::vector<std::string> out;
  for (size_t p : percenile_indices) {
    CHECK(p <= 100);
    out.emplace_back(std::to_string(percentiles[p]));
  }

  return nc::Join(out, " ");
}

static void ParseMatrix(const std::string& tm_id, const ctr::TrafficMatrix& tm,
                        size_t seed, ctr::Optimizer* optimizer,
                        ctr::RoutingConfigInfo* output_info,
                        ctr::RoutingConfigDeltaInfo* delta_info) {
  std::unique_ptr<ctr::RoutingConfiguration> baseline = optimizer->Optimize(tm);
  for (size_t i = 0; i < FLAGS_try_count; ++i) {
    std::mt19937 rnd(seed + i);
    auto rnd_tm = tm.Randomize(FLAGS_demand_fraction, FLAGS_flow_count_fraction,
                               FLAGS_aggregate_count, &rnd);
    if (!rnd_tm->ToDemandMatrix()->IsFeasible({})) {
      continue;
    }

    auto output = optimizer->Optimize(*rnd_tm);
    ctr::RoutingConfigurationDelta routing_config_delta =
        baseline->GetDifference(*output);
    delta_info->Add(routing_config_delta);
    output_info->Add(*output);
  }
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

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_topology_file.empty());

  std::vector<std::string> matrix_files;
  std::vector<std::string> to_glob = nc::Split(FLAGS_matrices, ",");
  for (const std::string string_to_glob : to_glob) {
    std::vector<std::string> globbed = nc::Glob(string_to_glob);
    matrix_files.insert(matrix_files.end(), globbed.begin(), globbed.end());
  }
  CHECK(!matrix_files.empty());

  std::vector<std::string> nodes_in_order;
  nc::net::GraphBuilder graph_builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology_file), &nodes_in_order);
  nc::net::GraphStorage graph(graph_builder);

  ctr::RoutingConfigInfo output_info;
  ctr::RoutingConfigDeltaInfo delta_info;

  for (const std::string& matrix_file : matrix_files) {
    auto demand_matrix = nc::lp::DemandMatrix::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(matrix_file), nodes_in_order, &graph);
    ctr::TrafficMatrix traffic_matrix(*demand_matrix,
                                      GetFlowCountMap(*demand_matrix));

    ctr::PathProvider path_provider(&graph);
    std::unique_ptr<ctr::Optimizer> optimizer;
    if (FLAGS_opt == "CTR") {
      optimizer = nc::make_unique<ctr::CTROptimizer>(&path_provider);
    } else if (FLAGS_opt == "MinMax") {
      optimizer = nc::make_unique<ctr::MinMaxOptimizer>(&path_provider);
    } else if (FLAGS_opt == "B4") {
      optimizer = nc::make_unique<ctr::B4Optimizer>(&path_provider);
    } else {
      LOG(FATAL) << "Unknown optimizer";
    }

    ParseMatrix(nc::StrCat(FLAGS_topology_file, ":", matrix_file),
                traffic_matrix, FLAGS_seed, optimizer.get(), &output_info,
                &delta_info);
  }

  nc::viz::NpyArray::Combine(output_info, delta_info).ToDisk("test_out");
}
