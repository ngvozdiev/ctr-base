#include <gflags/gflags.h>
#include <algorithm>
#include <chrono>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/algorithm.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/thread_runner.h"
#include "ncode_common/src/viz/grapher.h"
#include "demand_matrix_input.h"
#include "topology_input.h"

DEFINE_string(output_template, "tm_stat_template.rst",
              "RST template for the output");

using namespace std::chrono;

static void PlotCDF(const std::string& output, const std::string& x_label,
                    const std::string& title, const std::vector<double>& data) {
  nc::viz::DataSeries1D series;
  series.data = data;

  nc::viz::PlotParameters1D params;
  params.data_label = x_label;
  params.title = title;
  nc::viz::CDFPlot plot(params);
  plot.AddData(series);
  std::string svg = plot.PlotToSVG();
  LOG(INFO) << "S " << svg;
  nc::File::WriteStringToFileOrDie(svg, output);
}

static void PlotLine(const std::string& output, const std::string& x_label,
                     const std::string& y_label, const std::string& title,
                     const std::vector<std::pair<double, double>>& data) {
  nc::viz::DataSeries2D series;
  series.data = data;

  nc::viz::PlotParameters2D params;
  params.title = title;
  params.x_label = x_label;
  params.y_label = y_label;
  nc::viz::LinePlot plot(params);
  plot.AddData(series);
  plot.PlotToDir(output);
}

static void PlotDemandStats(const ctr::DemandMatrixAndFilename& input) {
  std::vector<double> demand_distances_ms;
  std::vector<double> demand_sizes_Mbps;
  std::vector<double> sp_utilizations;
  std::vector<std::pair<double, double>> cumulative_distances;

  const nc::lp::DemandMatrix& demand_matrix = *input.demand_matrix;
  const nc::net::GraphStorage* graph = demand_matrix.graph();
  nc::net::AllPairShortestPath sp({}, graph->AdjacencyList(), nullptr, nullptr);
  for (const auto& element : demand_matrix.elements()) {
    double distance_ms =
        duration_cast<milliseconds>(sp.GetDistance(element.src, element.dst))
            .count();
    double demand_Mbps = element.demand.Mbps();

    demand_distances_ms.emplace_back(distance_ms);
    demand_sizes_Mbps.emplace_back(demand_Mbps);
    cumulative_distances.emplace_back(distance_ms, demand_Mbps);
  }

  std::sort(cumulative_distances.begin(), cumulative_distances.end());
  double total = 0;
  for (auto& distance_and_demand : cumulative_distances) {
    double& demand = distance_and_demand.second;
    total += demand;
    demand = total;
  }

  nc::net::GraphLinkMap<double> utilizations = demand_matrix.SPUtilization();
  for (const auto& link_and_utilization : utilizations) {
    sp_utilizations.emplace_back(*link_and_utilization.second);
  }

  PlotCDF("demand_sizes.svg", "size (Mbps)", "CDF of demands\\' volumes",
          demand_sizes_Mbps);
  //
  //  std::string output_location =
  //      nc::Substitute(FLAGS_stats_pattern.c_str(), topology_filename,
  //                     FLAGS_locality, FLAGS_min_scale_factor, id);
  //  nc::File::RecursivelyCreateDir(output_location, 0777);
  //  PlotCDF(nc::StrCat(output_location, "/demand_distances"), "distance (ms)",
  //          "CDF of distances of demands\\' shortest paths",
  //          demand_distances_ms);
  //  PlotCDF(nc::StrCat(output_location, "/sp_utilizations"), "link
  //  utilization",
  //          "CDF of link utilizations when all aggregates are routed on the "
  //          "shortest path",
  //          sp_utilizations);
  //  PlotCDF(nc::StrCat(output_location, "/demand_sizes"), "size (Mbps)",
  //          "CDF of demands\\' volumes", demand_sizes_Mbps);
  //  PlotLine(nc::StrCat(output_location, "/cumulative_demands"),
  //           "SP distance (ms)", "Cumulative size (Mbps)",
  //           "Cumulative distance vs aggregate size", cumulative_distances);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> demands;
  std::tie(topologies, demands) = ctr::GetDemandMatrixInputs();

  nc::RunInParallel<ctr::DemandMatrixAndFilename>(
      demands, [](const ctr::DemandMatrixAndFilename& input) {
        PlotDemandStats(input);
      }, 1);
}
