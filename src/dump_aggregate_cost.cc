#include <gflags/gflags.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/map_util.h>
#include <ncode/net/net_common.h>
#include <ncode/strutil.h>
#include <ncode/viz/grapher.h>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "common.h"
#include "demand_matrix_input.h"
#include "opt/opt.h"
#include "opt/path_provider.h"
#include "topology_input.h"

DEFINE_double(sp_fraction, 1.2, "How far from the SP a path can be");
DEFINE_double(link_fraction_limit, 0.8,
              "At least this much of the SP's links can be routed around");
DEFINE_string(opt, "SP", "The optimizer to plot");
DEFINE_uint64(threads, 2, "How many threads to use");

// Will return the fraction of aggregates that do not fit.
static double GetDatapointForRC(const ctr::RoutingConfiguration& rc) {
  double total_count = rc.demands().size();
  return rc.OverloadedAggregates() / total_count;
}

static double GetDatapointForTopology(const nc::net::GraphStorage& graph) {
  return ctr::GetFractionOfPairsAboveLinkFraction(graph, FLAGS_sp_fraction,
                                                  FLAGS_link_fraction_limit);
}

static std::string GetFilename(const std::string& tm_file,
                               const std::string opt_string) {
  std::string tm_base = nc::StringReplace(tm_file, ".demands", "", true);
  return nc::StrCat(tm_base, "_", opt_string, ".rc");
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::tie(topologies, matrices) = ctr::GetDemandMatrixInputs(true);

  std::map<std::string, const ctr::TopologyAndFilename*> topologies_by_name;
  std::map<std::string, double> topology_datapoints;
  for (const auto& topology : topologies) {
    topologies_by_name[topology.file] = &topology;

    LOG(INFO) << "Getting datapoint for " << topology.file;
    double datapoint = GetDatapointForTopology(*topology.graph);
    topology_datapoints[topology.file] = datapoint;
  }

  std::map<std::string, std::vector<double>> rc_datapoints;
  for (const auto& matrix : matrices) {
    const ctr::TopologyAndFilename* topology_and_filename =
        nc::FindOrDieNoPrint(topologies_by_name, matrix.topology_file);

    ctr::PathProvider path_provider(topology_and_filename->graph.get());
    const nc::lp::DemandMatrix* demand_matrix = matrix.demand_matrix.get();
    std::string rc_filename = GetFilename(matrix.file, FLAGS_opt);
    if (!nc::File::Exists(rc_filename)) {
      LOG(INFO) << "Missing " << rc_filename;
      continue;
    }

    std::string rc_serialized = nc::File::ReadFileToStringOrDie(rc_filename);
    std::unique_ptr<ctr::TrafficMatrix> tm =
        ctr::TrafficMatrix::DistributeFromDemandMatrix(*demand_matrix);
    LOG(INFO) << "Will parse " << matrix.file << " at " << matrix.topology_file
              << " opt " << FLAGS_opt;
    auto rc = ctr::RoutingConfiguration::LoadFromSerializedText(
        *tm, topology_and_filename->node_order, rc_serialized, &path_provider);
    double datapoint = GetDatapointForRC(*rc);
    rc_datapoints[matrix.topology_file].emplace_back(datapoint);
  }

  std::vector<std::pair<double, double>> to_plot;
  for (const auto& topology_and_datapoints : rc_datapoints) {
    double topology_datapoint =
        nc::FindOrDie(topology_datapoints, topology_and_datapoints.first);
    for (double datapoint : topology_and_datapoints.second) {
      to_plot.emplace_back(topology_datapoint, datapoint);
    }
  }

  nc::viz::LinePlot plot;
  plot.TurnIntoScatterPlot();
  plot.AddData("", to_plot);
  plot.PlotToDir("scatter_out");
}
