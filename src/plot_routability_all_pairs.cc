#include <gflags/gflags.h>
#include <ncode/logging.h>
#include <ncode/viz/grapher.h>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "common.h"
#include "opt/opt.h"
#include "topology_input.h"

DEFINE_double(sp_fraction, 1.4, "How far from the SP a path can be");
DEFINE_string(output, "routability_all_pairs", "Directory to store the output");

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(FLAGS_output != "") << "Need output";

  std::vector<ctr::TopologyAndFilename> topologies = ctr::GetTopologyInputs();
  nc::viz::CDFPlot plot;
  for (const auto& topology : topologies) {
    std::vector<double> fractions;
    std::map<ctr::AggregateId, double> link_fraction_at_delay =
        ctr::GetLinkFractionAtDelay(*topology.graph, FLAGS_sp_fraction);
    for (const auto& aggregate_and_fraction : link_fraction_at_delay) {
      double fraction = aggregate_and_fraction.second;
      fractions.emplace_back(fraction);
    }
    plot.AddData("", fractions);
  }

  plot.PlotToDir(FLAGS_output);
}
