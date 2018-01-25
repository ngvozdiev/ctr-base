#include <gflags/gflags.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/map_util.h>
#include <ncode/net/net_common.h>
#include <ncode/strutil.h>
#include <ncode/viz/grapher.h>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "demand_matrix_input.h"
#include "opt/opt.h"
#include "plot_algorithm_eval_tools.h"
#include "topology_input.h"

DEFINE_string(
    optimizer_labels,
    "MinMax (10 paths),LDR (no flow counts),AB4,MinMax (low delay),LDR",
    "Labels of optimizers to plot.");
DEFINE_string(output_prefix, "sig",
              "Prefix that will be added to all output directories.");
DEFINE_double(
    routability_threshold, 0.5,
    "Threshold for routability, will do separate plots for topologies "
    "below and above that");
DEFINE_double(sp_fraction, 1.4, "How far from the SP a path can be");
DEFINE_double(link_fraction_limit, 0.7,
              "At least this much of the SP's links can be routed around");

using namespace ctr;
using namespace std::chrono;
using namespace ctr::alg_eval;

// Groups a series of TMSummary by optimizer. Will only consider graphs in the
// given set.
std::map<std::string, std::vector<const RCSummary*>> GroupByOpt(
    std::set<const nc::net::GraphStorage*>& graphs,
    const std::vector<std::unique_ptr<TMSummary>>& tm_summaries) {
  std::map<std::string, std::vector<const RCSummary*>> out;
  for (const auto& tm_summary : tm_summaries) {
    const nc::net::GraphStorage* graph = tm_summary->demand_matrix->graph();
    if (!nc::ContainsKey(graphs, graph)) {
      continue;
    }

    for (const auto& opt_and_summary : tm_summary->rcs) {
      out[opt_and_summary.first] = opt_and_summary.second.get();
    }
  }

  return out;
}

void PlotStretch(
    const std::string& prefix,
    const std::map<std::string, std::vector<const RCSummary*>>& summaries) {
  MultiRCSummaryPlotPack plot_pack;
  for (const auto& opt_and_summaries : summaries) {
    const std::string& opt = opt_and_summaries.first;
    const std::vector<const RCSummary*>& summaries_for_opt =
        opt_and_summaries.second;

    plot_pack.PlotStretchDistribution(opt, summaries_for_opt);
  }

  plot_pack.path_stretch_max_rel_plot().PlotToDir(
      nc::StrCat(prefix, "max_stretch"));
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::vector<std::unique_ptr<TMSummary>> summaries;
  std::tie(topologies, matrices, summaries) = GetTMSummaries();

  // Need to know the routability paths for each graph.
  std::set<const nc::net::GraphStorage*> high_routability_graphs;
  std::set<const nc::net::GraphStorage*> low_routability_graphs;
  for (const auto& topology : topologies) {
    const nc::net::GraphStorage* graph = topology.graph.get();
    double routability = GetFractionOfPairsAboveLinkFraction(
        *graph, FLAGS_sp_fraction, FLAGS_link_fraction_limit);

    if (routability > FLAGS_routability_threshold) {
      high_routability_graphs.insert(graph);
    } else {
      low_routability_graphs.insert(graph);
    }
  }

  auto hr_by_opt = GroupByOpt(high_routability_graphs, summaries);
  auto lr_by_opt = GroupByOpt(low_routability_graphs, summaries);

  std::string hr_prefix = nc::StrCat(FLAGS_output_prefix, "_hr_");
  std::string lr_prefix = nc::StrCat(FLAGS_output_prefix, "_lr_");
  PlotStretch(hr_prefix, hr_by_opt);
  PlotStretch(lr_prefix, lr_by_opt);
}