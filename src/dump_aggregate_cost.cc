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
DEFINE_uint64(threads, 2, "How many threads to use");

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::viz::CDFPlot fraction_plot;
  std::vector<ctr::TopologyAndFilename> topologies = ctr::GetTopologyInputs();
  for (const auto& topology_and_filename : topologies) {
    const nc::net::GraphStorage& graph = *topology_and_filename.graph;
    std::map<ctr::AggregateId, double> fractions_at_delay =
        ctr::GetLinkFractionAtDelay(graph, FLAGS_sp_fraction);
    std::vector<double> fractions;
    for (const auto& aggregate_and_fraction : fractions_at_delay) {
      fractions.emplace_back(aggregate_and_fraction.second);
    }

    std::string label = nc::File::ExtractFileName(topology_and_filename.file);
    label = nc::StrCat(label, " n", static_cast<uint32_t>(graph.NodeCount()),
                       " e", static_cast<uint32_t>(graph.LinkCount()));
    fraction_plot.AddData(label, fractions);
  }

  fraction_plot.PlotToDir("fractions_out");
}
