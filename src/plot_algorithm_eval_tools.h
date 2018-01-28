// Classes and utilities to parse the metrics generated by opt_eval_util.

#include <ncode/common.h>
#include <ncode/net/net_common.h>
#include <ncode/viz/grapher.h>
#include <iostream>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

#include "demand_matrix_input.h"
#include "topology_input.h"

namespace ctr {
namespace alg_eval {

// Graph, seed and load and locality.
using TMKey = std::tuple<const nc::net::GraphStorage*, double, double, double>;

struct TMSummary;

struct RCSummary {
  RCSummary() {}

  double total_sp_delay = 0;
  double total_delay = 0;

  // One value per path.
  std::vector<double> abs_stretches;
  std::vector<double> rel_stretches;
  std::vector<double> flow_counts;
  std::vector<bool> overloaded;

  // One value per aggregate.
  std::vector<double> path_counts;

  // Link utilizations.
  nc::net::GraphLinkMap<double> utilizations;

  // The optimizer that was used to generate the summary.
  const std::string* opt = nullptr;

  // The parent.
  const TMSummary* parent = nullptr;

  // How long optimization took.
  double runtime_ms = 0;

  DISALLOW_COPY_AND_ASSIGN(RCSummary);
};

struct TMSummary {
  TMSummary() {}

  // The demand matrix.
  const nc::lp::DemandMatrix* demand_matrix = nullptr;

  // Identifies the graph, seed, load and locality.
  TMKey key;

  // File names.
  const std::string* topology_file_name = nullptr;
  const std::string* demand_matrix_file_name = nullptr;

  // RCs for this traffic matrix.
  std::map<std::string, std::unique_ptr<RCSummary>> rcs;

  DISALLOW_COPY_AND_ASSIGN(TMSummary);
};

class InterestingTable {
 public:
  InterestingTable(const std::vector<std::string>& header) : header_(header) {}

  void Add(const std::vector<std::string>& row, const TMSummary* tm);

  void Add(const std::vector<std::string>& row);

  const std::vector<const TMSummary*>& tms() const { return tms_; }

  const std::vector<std::vector<std::string>>& rows() const { return rows_; }

  std::string ToRSTTable(
      std::vector<const TMSummary*>* interestring_tm_list) const;

 private:
  std::string title_;
  std::vector<std::string> header_;
  std::vector<std::vector<std::string>> rows_;
  std::vector<const TMSummary*> tms_;
};

class MultiRCSummaryPlotPack {
 public:
  MultiRCSummaryPlotPack();

  void PlotStretchDistribution(const std::string& optimizer,
                               const std::vector<const RCSummary*>& rcs);

  void PlotPathRatios(const std::string& optimizer,
                      const std::vector<const RCSummary*>& rcs);

  void PlotPathCounts(const std::string& optimizer,
                      const std::vector<const RCSummary*>& rcs);

  void PlotLinkUtilizations(const std::string& optimizer,
                            const std::vector<const RCSummary*>& rcs);

  void PlotLinkScalesAtDelay(const std::vector<const TMSummary*>& all_demands);

  const InterestingTable& aggregate_path_count_interesting() const {
    return aggregate_path_count_interesting_;
  }

  const nc::viz::CDFPlot& aggregate_path_count_plot() const {
    return aggregate_path_count_plot_;
  }

  const nc::viz::CDFPlot& link_scales_plot() const { return link_scales_plot_; }

  const InterestingTable& link_utilization_interesting() const {
    return link_utilization_interesting_;
  }

  const nc::viz::CDFPlot& link_utilization_plot() const {
    return link_utilization_plot_;
  }

  const InterestingTable& path_stretch_max_rel_interesting() const {
    return path_stretch_max_rel_interesting_;
  }

  const nc::viz::CDFPlot& path_stretch_max_rel_plot() const {
    return path_stretch_max_rel_plot_;
  }

  const nc::viz::LinePlot& path_stretch_rel_plot() const {
    return path_stretch_rel_plot_;
  }

  const InterestingTable& ratios_interesting() const {
    return ratios_interesting_;
  }

  const InterestingTable& link_scales_interesting() const {
    return link_scales_interesting_;
  }

  const nc::viz::CDFPlot& ratios_plot() const { return ratios_plot_; }

 private:
  nc::viz::LinePlot path_stretch_rel_plot_;

  nc::viz::CDFPlot path_stretch_max_rel_plot_;
  InterestingTable path_stretch_max_rel_interesting_;

  nc::viz::CDFPlot link_utilization_plot_;
  InterestingTable link_utilization_interesting_;

  nc::viz::CDFPlot ratios_plot_;
  InterestingTable ratios_interesting_;

  nc::viz::CDFPlot aggregate_path_count_plot_;
  InterestingTable aggregate_path_count_interesting_;

  nc::viz::CDFPlot link_scales_plot_;
  InterestingTable link_scales_interesting_;
};

class SingleRCSummaryPlotPack {
 public:
  SingleRCSummaryPlotPack();

  void PlotCumulativeDistances(const nc::lp::DemandMatrix& demand_matrix,
                               const nc::net::AllPairShortestPath& sp);

  void PlotDemandSizes(const nc::lp::DemandMatrix& demand_matrix);

  void PlotSPUtilizations(const nc::lp::DemandMatrix& demand_matrix);

  void PlotTotalDelayAtLinkScale(const nc::lp::DemandMatrix& demand_matrix);

  void PlotDelayVsUtilization(const std::string& opt_label,
                              const RCSummary& rc);

  void PlotAbsoluteStretch(const std::string& opt_label, const RCSummary& rc);

  const nc::viz::LinePlot& cumulative_demands_plot() const {
    return cumulative_demands_plot_;
  }

  const nc::viz::LinePlot& cumulative_demands_hop_plot() const {
    return cumulative_demands_hop_plot_;
  }

  const InterestingTable& demand_sizes_interesting() const {
    return demand_sizes_interesting_;
  }

  const nc::viz::CDFPlot& demand_sizes_plot() const {
    return demand_sizes_plot_;
  }

  const InterestingTable& sp_utilizations_interesting() const {
    return sp_utilizations_interesting_;
  }

  const nc::viz::CDFPlot& sp_utilizations_plot() const {
    return sp_utilizations_plot_;
  }

  const nc::viz::LinePlot& total_delay_at_link_scale_plot() const {
    return total_delay_at_link_scale_plot_;
  }

  const nc::viz::LinePlot& link_delay_vs_link_utilization_plot() const {
    return link_delay_vs_link_utilization_plot_;
  }

  const nc::viz::LinePlot& absolute_path_stretch_plot() const {
    return absolute_path_stretch_plot_;
  }

 private:
  // Plots demand sizes for the topology.
  nc::viz::CDFPlot demand_sizes_plot_;
  InterestingTable demand_sizes_interesting_;

  // Plots link utilizations when everything is routed on the SP.
  nc::viz::CDFPlot sp_utilizations_plot_;
  InterestingTable sp_utilizations_interesting_;

  // Plots cumulative distance.
  nc::viz::LinePlot cumulative_demands_plot_;
  nc::viz::LinePlot cumulative_demands_hop_plot_;

  // Plots how adding headroom to all links affects the total delay experienced
  // by all flows in the topology.
  nc::viz::LinePlot total_delay_at_link_scale_plot_;

  // Plot with one curve per optimizer. X axis is links, Y axis is utilization.
  nc::viz::LinePlot link_delay_vs_link_utilization_plot_;

  // Plot with one curve per optimizer. X axis is per-flow path stretch in
  // milliseconds.
  nc::viz::LinePlot absolute_path_stretch_plot_;
};

std::tuple<std::vector<TopologyAndFilename>,
           std::vector<DemandMatrixAndFilename>,
           std::vector<std::unique_ptr<TMSummary>>>
GetTMSummaries();

}  // namespace alg_eval
}  // namespace ctr
