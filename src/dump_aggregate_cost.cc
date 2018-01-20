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
DEFINE_string(opt, "SP,B4,MinMaxK10,CTR", "The optimizers to plot");
DEFINE_bool(plot_delay, true, "If true will plot delay, if false capacity");

// Will return the fraction of aggregates that do not fit.
static double GetDatapointForRC(const ctr::RoutingConfiguration& rc) {
  if (FLAGS_plot_delay) {
    nc::net::Delay total = rc.TotalPerFlowDelay(false);
    nc::net::Delay total_sp = rc.TotalPerFlowDelay(true);
    double fraction = static_cast<double>(total.count()) / total_sp.count();
    return fraction;
  }

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

static void ParseOpt(
    const std::string& opt,
    const std::vector<ctr::DemandMatrixAndFilename>& matrices,
    const std::map<std::string, const ctr::TopologyAndFilename*>&
        topologies_by_name,
    const std::map<std::string, double>& topology_datapoints) {
  std::map<std::string, std::vector<double>> rc_datapoints;
  for (const auto& matrix : matrices) {
    const ctr::TopologyAndFilename* topology_and_filename =
        nc::FindOrDieNoPrint(topologies_by_name, matrix.topology_file);

    ctr::PathProvider path_provider(topology_and_filename->graph.get());
    const nc::lp::DemandMatrix* demand_matrix = matrix.demand_matrix.get();
    std::string rc_filename = GetFilename(matrix.file, opt);
    if (!nc::File::Exists(rc_filename)) {
      LOG(INFO) << "Missing " << rc_filename;
      continue;
    }

    std::string rc_serialized = nc::File::ReadFileToStringOrDie(rc_filename);
    std::unique_ptr<ctr::TrafficMatrix> tm =
        ctr::TrafficMatrix::DistributeFromDemandMatrix(*demand_matrix);
    LOG(INFO) << "Will parse " << matrix.file << " at " << matrix.topology_file
              << " opt " << opt;
    auto rc = ctr::RoutingConfiguration::LoadFromSerializedText(
        *tm, topology_and_filename->node_order, rc_serialized, &path_provider);
    double datapoint = GetDatapointForRC(*rc);
    rc_datapoints[matrix.topology_file].emplace_back(datapoint);
  }

  // Will only label the top 10.
  std::vector<std::pair<double, std::string>> topology_values_vector;
  for (const auto& topology_and_datapoints : rc_datapoints) {
    const std::string& topology = topology_and_datapoints.first;
    double value = *std::max_element(topology_and_datapoints.second.begin(),
                                     topology_and_datapoints.second.end());
    topology_values_vector.emplace_back(value, topology);
  }

  std::sort(topology_values_vector.begin(), topology_values_vector.end(),
            std::greater<std::pair<double, std::string>>());
  topology_values_vector.resize(std::min(topology_values_vector.size(), 7ul));

  std::set<std::string> topologies_to_label;
  for (const auto& value_and_topology : topology_values_vector) {
    topologies_to_label.emplace(value_and_topology.second);
  }

  std::string x_label = nc::Substitute(
      "Fraction of pairs with alternative path for at least $0% of SP links, "
      "$1% of SP delay",
      static_cast<int>(FLAGS_link_fraction_limit * 100),
      static_cast<int>(FLAGS_sp_fraction * 100));

  std::string y_label = FLAGS_plot_delay
                            ? "total float delay / total sp flow delay"
                            : "fraction of pairs that cross overloaded links";
  nc::viz::LinePlot plot({"", x_label, y_label});

  plot.TurnIntoScatterPlot();
  std::vector<std::pair<double, double>> rest;
  std::vector<std::pair<double, double>> medians;
  std::vector<std::pair<double, double>> means;
  for (auto& topology_and_datapoints : rc_datapoints) {
    double topology_datapoint =
        nc::FindOrDie(topology_datapoints, topology_and_datapoints.first);
    const std::string& topology_name = topology_and_datapoints.first;
    std::string name_stripped = nc::File::ExtractFileName(topology_name);

    std::vector<double>& datapoints = topology_and_datapoints.second;
    std::sort(datapoints.begin(), datapoints.end());
    double med = datapoints[datapoints.size() / 2];
    medians.emplace_back(topology_datapoint, med);

    double mean = std::accumulate(datapoints.begin(), datapoints.end(), 0.0) /
                  datapoints.size();
    means.emplace_back(topology_datapoint, mean);

    std::vector<std::pair<double, double>> to_plot;
    for (double datapoint : datapoints) {
      to_plot.emplace_back(topology_datapoint, datapoint);
    }

    if (nc::ContainsKey(topologies_to_label, topology_name)) {
      plot.AddData(name_stripped, to_plot);
    } else {
      rest.insert(rest.end(), to_plot.begin(), to_plot.end());
    }
  }

  plot.AddData("Rest", rest);
  plot.AddData("Medians", medians);
  plot.AddData("Means", means);

  std::string suffix = FLAGS_plot_delay ? "delay" : "capacity";
  plot.PlotToDir(nc::StrCat("scatter_out_", opt, "_", suffix));
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::tie(topologies, matrices) = ctr::GetDemandMatrixInputs(false);

  std::map<std::string, const ctr::TopologyAndFilename*> topologies_by_name;
  std::map<std::string, double> topology_datapoints;
  for (const auto& topology : topologies) {
    topologies_by_name[topology.file] = &topology;

    LOG(INFO) << "Getting datapoint for " << topology.file;
    double datapoint = GetDatapointForTopology(*topology.graph);
    topology_datapoints[topology.file] = datapoint;
  }

  std::vector<std::string> opts = nc::Split(FLAGS_opt, ",", true);
  for (const std::string& opt : opts) {
    ParseOpt(opt, matrices, topologies_by_name, topology_datapoints);
  }
}
