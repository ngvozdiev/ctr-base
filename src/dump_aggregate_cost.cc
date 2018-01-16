#include <gflags/gflags.h>
#include <ncode/logging.h>
#include <ncode/net/algorithm.h>
#include <ncode/net/net_common.h>
#include <ncode/viz/grapher.h>
#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <map>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "common.h"
#include "opt/opt.h"
#include "topology_input.h"

DEFINE_uint64(seed, 1ul, "Seed for the generated TM.");
DEFINE_uint64(k, 10000ul, "Number of paths");

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::vector<ctr::TopologyAndFilename> topologies = ctr::GetTopologyInputs();
  CHECK(topologies.size() == 1) << "Too many topologies";

  const nc::net::GraphStorage* graph = topologies[0].graph.get();
  std::vector<nc::net::GraphNodeIndex> nodes_ordered;
  for (nc::net::GraphNodeIndex node : graph->AllNodes()) {
    nodes_ordered.emplace_back(node);
  }

  std::map<ctr::AggregateId, double> capacity_at_delay =
      ctr::GetCapacityAtDelay(*graph, 1.2);
  std::vector<double> capacities;
  for (const auto& aggregate_and_capacity : capacity_at_delay) {
    capacities.emplace_back(aggregate_and_capacity.second);
  }
  nc::viz::CDFPlot capacities_plot;
  capacities_plot.AddData("", capacities);
  capacities_plot.PlotToDir("capacities_out");

  std::mt19937 rnd(FLAGS_seed);
  std::shuffle(nodes_ordered.begin(), nodes_ordered.end(), rnd);
  CHECK(nodes_ordered.size() > 1);

  nc::net::GraphNodeIndex src = nodes_ordered[0];
  nc::net::GraphNodeIndex dst = nodes_ordered[1];

  std::vector<std::pair<double, double>> delays;
  nc::net::KShortestPathsGenerator ksp_gen(src, dst, *graph, {});
  for (size_t i = 0; i < FLAGS_k; ++i) {
    const nc::net::Walk* p = ksp_gen.KthShortestPathOrNull(i);
    if (p == nullptr) {
      break;
    }

    double cost = p->delay().count();
    delays.emplace_back(i, cost);
  }

  nc::viz::LinePlot plot;
  plot.AddData("original", delays);
  plot.PlotToDir("dump_aggregate_out");
}
