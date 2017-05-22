#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/strutil.h"
#include "../common.h"
#include "ctr.h"
#include "path_provider.h"

DEFINE_string(topology_file, "", "A topology file");
DEFINE_string(matrices, "", "One or more matrices");
DEFINE_double(demand_fraction, 0.05, "By how much to vary demand");
DEFINE_double(flow_count_fraction, 0, "By how much to vary flow counts");
DEFINE_uint64(aggregate_count, 1, "How many aggregates to change");
DEFINE_uint64(try_count, 1, "How many tries to perform");

static double MaxFractionDelta(
    const std::map<ctr::AggregateId, ctr::AggregateDelta>& difference_map) {
  double max_fraction_delta = 0;
  for (const auto& aggregate_and_delta : difference_map) {
    const ctr::AggregateDelta& delta = aggregate_and_delta.second;
    max_fraction_delta = std::max(max_fraction_delta, delta.fraction_delta);
  }

  return max_fraction_delta;
}

static void ParseMatrix(const ctr::TrafficMatrix& tm, size_t seed,
                        ctr::Optimizer* optimizer) {
  std::unique_ptr<ctr::RoutingConfiguration> baseline = optimizer->Optimize(tm);
  ctr::RoutingConfigurationDelta worst_try;
  double worst_delta = 0;

  std::mt19937 rnd(seed);
  for (size_t i = 0; i < FLAGS_try_count; ++i) {
    auto rnd_tm = tm.Randomize(FLAGS_demand_fraction, FLAGS_flow_count_fraction,
                               FLAGS_aggregate_count, &rnd);
    auto output = optimizer->Optimize(*rnd_tm);
    ctr::RoutingConfigurationDelta routing_config_delta =
        baseline->GetDifference(*output);

    double delta = routing_config_delta.total_flow_fraction_delta;
    if (delta > worst_delta) {
      worst_delta = delta;
      worst_try = routing_config_delta;
    }
  }

  LOG(ERROR) << "Worst try delta " << worst_delta << " mf "
             << MaxFractionDelta(worst_try.aggregates);
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

  for (const std::string& matrix_file : matrix_files) {
    auto demand_matrix = nc::lp::DemandMatrix::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(matrix_file), nodes_in_order, &graph);
    ctr::TrafficMatrix traffic_matrix(*demand_matrix,
                                      GetFlowCountMap(*demand_matrix));

    auto path_provider = nc::make_unique<ctr::PathProvider>(&graph);
    ctr::CTROptimizer optimizer(std::move(path_provider));
    ParseMatrix(traffic_matrix, 1, &optimizer);
  }
}
