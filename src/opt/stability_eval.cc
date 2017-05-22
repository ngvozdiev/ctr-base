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
DEFINE_string(opt, "CTR", "Optimizer to run. One of CTR,MinMax,B4");

static void ParseMatrix(const std::string& tm_id, const ctr::TrafficMatrix& tm,
                        size_t seed, ctr::Optimizer* optimizer) {
  std::unique_ptr<ctr::RoutingConfiguration> baseline = optimizer->Optimize(tm);
  std::vector<ctr::AggregateDelta> deltas;
  double worst_total_demand_delta = 0;
  double worst_total_flow_delta = 0;
  std::vector<double> fraction_deltas;

  std::mt19937 rnd(seed);
  for (size_t i = 0; i < FLAGS_try_count; ++i) {
    auto rnd_tm = tm.Randomize(FLAGS_demand_fraction, FLAGS_flow_count_fraction,
                               FLAGS_aggregate_count, &rnd);
    if (!rnd_tm->ToDemandMatrix()->IsFeasible({})) {
      continue;
    }

    auto output = optimizer->Optimize(*rnd_tm);
    ctr::RoutingConfigurationDelta routing_config_delta =
        baseline->GetDifference(*output);

    double demand_delta = routing_config_delta.total_volume_fraction_delta;
    double flow_delta = routing_config_delta.total_flow_fraction_delta;
    worst_total_demand_delta = std::max(worst_total_demand_delta, demand_delta);
    worst_total_flow_delta = std::max(worst_total_flow_delta, flow_delta);

    for (const auto& aggregate_id_and_delta : routing_config_delta.aggregates) {
      const ctr::AggregateDelta& aggregate_delta =
          aggregate_id_and_delta.second;
      fraction_deltas.emplace_back(aggregate_delta.fraction_delta);
    }
  }

  std::vector<double> delta_percentiles = nc::Percentiles(&fraction_deltas);
  std::cout << tm_id << " " << worst_total_demand_delta << " "
            << tm.demands().size() << " " << worst_total_flow_delta << " "
            << delta_percentiles[50] << " " << delta_percentiles[90] << " "
            << delta_percentiles[95] << " " << delta_percentiles[100] << "\n";
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
    std::unique_ptr<ctr::Optimizer> optimizer;
    if (FLAGS_opt == "CTR") {
      optimizer = nc::make_unique<ctr::CTROptimizer>(std::move(path_provider));
    } else if (FLAGS_opt == "MinMax") {
      optimizer =
          nc::make_unique<ctr::MinMaxOptimizer>(std::move(path_provider));
    } else if (FLAGS_opt == "B4") {
      optimizer = nc::make_unique<ctr::B4Optimizer>(std::move(path_provider));
    } else {
      LOG(FATAL) << "Unknown optimizer";
    }

    ParseMatrix(nc::StrCat(FLAGS_topology_file, ":", matrix_file),
                traffic_matrix, 1, optimizer.get());
  }
}
