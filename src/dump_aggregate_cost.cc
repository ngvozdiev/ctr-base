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

DEFINE_string(sp_fractions, "1.2", "How far from the SP a path can be");

static std::vector<double> GetCommaSeparated(const std::string& value) {
  std::vector<std::string> split = nc::Split(value, ",");
  std::vector<double> out;
  for (const auto& piece : split) {
    double v;
    CHECK(nc::safe_strtod(piece, &v));
    out.emplace_back(v);
  }

  return out;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::vector<ctr::TopologyAndFilename> topologies = ctr::GetTopologyInputs();
  CHECK(topologies.size() == 1) << "Too many topologies";

  nc::viz::CDFPlot capacities_plot;
  const nc::net::GraphStorage* graph = topologies[0].graph.get();
  std::vector<double> values = GetCommaSeparated(FLAGS_sp_fractions);
  for (double v : values) {
    std::map<ctr::AggregateId, uint64_t> counts_at_delay =
        ctr::GetPathCountsAtDelay(*graph, v);
    std::vector<double> capacities;
    for (const auto& aggregate_and_capacity : counts_at_delay) {
      capacities.emplace_back(aggregate_and_capacity.second);
    }

    capacities_plot.AddData(nc::StrCat("fraction ", v), capacities);
  }

  capacities_plot.PlotToDir("capacities_out");
}
