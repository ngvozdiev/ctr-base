#include <gflags/gflags.h>
#include <ncode/file.h>
#include <ncode/viz/grapher.h>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "common.h"
#include "topology_input.h"
#include "opt/opt.h"

DEFINE_double(sp_fraction, 1.2, "How far from the SP a path can be");

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::viz::CDFPlot path_counts_plot;
  std::vector<ctr::TopologyAndFilename> topologies = ctr::GetTopologyInputs();
  for (const auto& topology_and_filename : topologies) {
    const nc::net::GraphStorage& graph = *topology_and_filename.graph;
    std::map<ctr::AggregateId, uint64_t> counts_at_delay =
        ctr::GetPathCountsAtDelay(graph, FLAGS_sp_fraction);
    std::vector<double> path_counts;
    for (const auto& aggregate_and_path_count : counts_at_delay) {
      path_counts.emplace_back(aggregate_and_path_count.second);
    }

    std::string label = nc::File::ExtractFileName(topology_and_filename.file);
    path_counts_plot.AddData(label, path_counts);
  }

  path_counts_plot.PlotToDir("path_counts_out");
}
