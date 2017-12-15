#include <gflags/gflags.h>
#include <algorithm>
#include <chrono>
#include <iostream>
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
#include "ncode_common/src/viz/ctemplate/template.h"
#include "ncode_common/src/viz/ctemplate/template_dictionary.h"
#include "ncode_common/src/viz/ctemplate/template_enums.h"
#include "ncode_common/src/viz/grapher.h"
#include "demand_matrix_input.h"
#include "topology_input.h"

DEFINE_string(output_root, "./", "Where to save plots");

using namespace std::chrono;

static nc::viz::CDFPlot GetCDFPlot(const std::string& x_label,
                                   const std::string& title) {
  nc::viz::PlotParameters1D params;
  params.data_label = x_label;
  params.title = title;
  nc::viz::CDFPlot plot(params);
  return plot;
}

static nc::viz::LinePlot GetLinePlot(const std::string& x_label,
                                     const std::string& y_label,
                                     const std::string& title) {
  nc::viz::PlotParameters2D params;
  params.title = title;
  params.x_label = x_label;
  params.y_label = y_label;
  nc::viz::LinePlot plot(params);
  return plot;
}

static nc::viz::HeatmapPlot GetHeatmapPlot(const std::string& x_label,
                                           const std::string& y_label,
                                           const std::string& title) {
  nc::viz::PlotParameters2D params;
  params.title = title;
  params.x_label = x_label;
  params.y_label = y_label;
  nc::viz::HeatmapPlot plot(params);
  return plot;
}

static void PlotDemandStats(
    const std::vector<ctr::DemandMatrixAndFilename>& inputs) {
  nc::viz::CDFPlot demand_sizes_plot = GetCDFPlot(
      "size (Mbps)", "Demand sizes (all possible demands considered)");
  nc::viz::CDFPlot demand_distances_plot =
      GetCDFPlot("SP length (ms)", "Length of each demand\\'s shortest path");
  nc::viz::CDFPlot sp_utilizations_plot =
      GetCDFPlot("link utilization", "Shortest path link utilization");
  nc::viz::LinePlot cumulative_demands_plot =
      GetLinePlot("SP distance (ms)", "Cumulative size (Mbps)",
                  "Cumulative distance vs aggregate size");

  ctemplate::TemplateDictionary dict("plot");
  for (size_t i = 0; i < inputs.size(); ++i) {
    const ctr::DemandMatrixAndFilename& input = inputs[i];

    std::vector<double> demand_distances_ms;
    std::vector<double> demand_sizes_Mbps;
    std::vector<double> sp_utilizations;
    std::vector<std::pair<double, double>> cumulative_distances;

    const nc::lp::DemandMatrix& demand_matrix = *input.demand_matrix;
    const nc::net::GraphStorage* graph = demand_matrix.graph();
    nc::net::AllPairShortestPath sp({}, graph->AdjacencyList(), nullptr,
                                    nullptr);
    for (const auto& element : demand_matrix.elements()) {
      double distance_ms =
          duration_cast<milliseconds>(sp.GetDistance(element.src, element.dst))
              .count();
      double demand_Mbps = element.demand.Mbps();
      if (demand_Mbps == 0) {
        continue;
      }

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

    size_t node_count = graph->NodeCount();
    double possible_count = node_count * (node_count - 1);
    CHECK(demand_sizes_Mbps.size() < possible_count);

    std::string locality =
        nc::FindOrDie(demand_matrix.properties(), "locality");
    std::string id = nc::StrCat("locality ", locality);

    demand_sizes_plot.AddData(id, demand_sizes_Mbps);
    demand_distances_plot.AddData(id, demand_distances_ms);
    sp_utilizations_plot.AddData(id, sp_utilizations);
    cumulative_demands_plot.AddData(id, cumulative_distances);

    ctemplate::TemplateDictionary* subdict =
        dict.AddSectionDictionary("summary_table_rows");

    size_t element_count = demand_matrix.elements().size();
    double mcsf = demand_matrix.MaxCommodityScaleFactor({}, 1.0);
    double element_count_fraction = element_count / possible_count;
    subdict->SetValue("name", id);
    subdict->SetValue("load_value", nc::ToStringMaxDecimals(mcsf, 2));
    subdict->SetValue("demand_count_value", std::to_string(element_count));
    subdict->SetValue("demand_fraction_value",
                      nc::ToStringMaxDecimals(element_count_fraction, 3));

    nc::viz::HeatmapPlot heatmap_plot = GetHeatmapPlot(
        "destination", "source", nc::StrCat("Demand sizes (", id, ")"));
    heatmap_plot.set_symlog(true);
    nc::net::GraphNodeMap<nc::net::GraphNodeMap<double>> values;

    for (const auto& element : demand_matrix.elements()) {
      double demand_Mbps = element.demand.Mbps();
      values[element.src][element.dst] = demand_Mbps;
    }

    // Will sort nodes based on outgoing demand.
    std::vector<std::pair<double, nc::net::GraphNodeIndex>> all_nodes;
    for (nc::net::GraphNodeIndex src_node : graph->AllNodes()) {
      nc::net::GraphNodeMap<double>& inner_values = values[src_node];
      double total = 0;
      for (nc::net::GraphNodeIndex dst_node : graph->AllNodes()) {
        double value = inner_values[dst_node];
        total += value;
      }

      all_nodes.emplace_back(total, src_node);
    }
    std::sort(all_nodes.begin(), all_nodes.end());

    for (const auto& demand_and_src_node : all_nodes) {
      nc::net::GraphNodeIndex src_node = demand_and_src_node.second;

      std::vector<double> row;
      nc::net::GraphNodeMap<double>& inner_values = values[src_node];
      for (const auto& demand_and_dst_node : all_nodes) {
        nc::net::GraphNodeIndex dst_node = demand_and_dst_node.second;
        double value = inner_values[dst_node];
        row.emplace_back(value);
      }

      heatmap_plot.AddData(row);
    }

    std::string plot_location = nc::StrCat(id, "_heatmap.svg");
    plot_location = nc::StringReplace(plot_location, " ", "_", true);

    heatmap_plot.PlotToSVGFile(plot_location);
    ctemplate::TemplateDictionary* heatmap_subdict =
        dict.AddSectionDictionary("demand_heatmaps");
    heatmap_subdict->SetValue("heatmap_location", plot_location);
  }

  demand_sizes_plot.PlotToSVGFile("demand_sizes.svg");
  dict.SetValue("demands_cdf_location", "demand_sizes.svg");

  demand_distances_plot.PlotToSVGFile("demand_distances.svg");
  dict.SetValue("delays_cdf_location", "demand_distances.svg");

  sp_utilizations_plot.PlotToSVGFile("sp_utilizations.svg");
  dict.SetValue("sp_utilization_location", "sp_utilizations.svg");

  cumulative_demands_plot.PlotToSVGFile("cumulative_demands.svg");
  dict.SetValue("cumulative_demands_location", "cumulative_demands.svg");
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> demands;
  std::tie(topologies, demands) = ctr::GetDemandMatrixInputs();

  PlotDemandStats(demands);
}
