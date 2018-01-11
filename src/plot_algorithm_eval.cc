#include <gflags/gflags.h>
#include <ncode/common.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/map_util.h>
#include <ncode/net/algorithm.h>
#include <ncode/net/net_common.h>
#include <ncode/stats.h>
#include <ncode/strutil.h>
#include <ncode/viz/ctemplate/template.h>
#include <ncode/viz/ctemplate/template_dictionary.h>
#include <ncode/viz/ctemplate/template_enums.h>
#include <ncode/viz/grapher.h>
#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "common.h"
#include "demand_matrix_input.h"
#include "opt/path_provider.h"
#include "topology_input.h"

DEFINE_string(optimizers, "MinMaxLD,CTRNFC,B4,CTR,MinMaxK10",
              "Optimizers to plot.");
DEFINE_string(output_prefix, "",
              "Prefix that will be added to all output directories.");

static constexpr char kTopLevelTemplate[] = "../rst/top_level_template.rst";
static constexpr char kSingleTopTemplate[] = "../rst/single_top_template.rst";
static constexpr char kSummaryTemplate[] = "../rst/summary_template.rst";
static constexpr char kLocalityLevelTemplate[] = "../rst/tm_stat_template.rst";
static constexpr char kSingleTmTemplate[] = "single_tm_template.rst";

static constexpr size_t kDiscreteMultiplier = 100;
static constexpr size_t kPercentilesCount = 10000;

using namespace ctr;
using namespace std::chrono;

using TopologyAndTM = std::pair<std::string, std::string>;

static double ToMs(nc::net::Delay delay) {
  return duration_cast<microseconds>(delay).count() / 1000.0;
}

// Graph, seed and load and locality.
using TMKey = std::tuple<const nc::net::GraphStorage*, double, double, double>;

struct RCSummary {
  double total_sp_delay = 0;
  double total_delay = 0;

  // One value per path.
  std::vector<double> abs_stretches;
  std::vector<double> rel_stretches;
  std::vector<double> flow_counts;

  // One value per aggregate.
  std::vector<double> path_counts;

  // Link utilizations for links that are less than the median link lenght and
  // more than it.
  std::vector<double> short_link_utilizations;
  std::vector<double> long_link_utilizations;

  // The demand matrix.
  const nc::lp::DemandMatrix* demand_matrix;

  // Identifies the graph, seed, load and locality.
  TMKey key;

  // The optimizer that was used to generate the summary.
  std::string opt;

  // File names.
  const std::string* topology_file_name;
  const std::string* demand_matrix_file_name;

  DISALLOW_COPY_AND_ASSIGN(RCSummary);
};

class InterestingTable {
 public:
  void Add(const std::string& field_one, double value, const RCSummary* rc) {
    std::vector<std::string> new_row;
    new_row.emplace_back(field_one);
    new_row.emplace_back(nc::StrCat(value));
    rows_.emplace_back(new_row);
    rcs_.emplace_back(rc);
  }

  void Add(const std::vector<std::string>& row) { rows_.emplace_back(row); }

  const std::vector<const RCSummary*>& rcs() const { return rcs_; }

  const std::vector<std::vector<std::string>>& rows() const { return rows_; }

 private:
  std::vector<std::vector<std::string>> rows_;
  std::vector<const RCSummary*> rcs_;
};

static nc::net::Delay MedianLinkDelay(const nc::net::GraphStorage& graph) {
  std::vector<nc::net::Delay> delays;
  for (nc::net::GraphLinkIndex link : graph.AllLinks()) {
    nc::net::Delay delay = graph.GetLink(link)->delay();
    delays.emplace_back(delay);
  }

  std::sort(delays.begin(), delays.end());
  CHECK(!delays.empty());
  return delays[delays.size() / 2];
}

static std::unique_ptr<RCSummary> ParseRC(
    const RoutingConfiguration& rc, const nc::net::AllPairShortestPath& sp,
    const nc::lp::DemandMatrix* demand_matrix) {
  auto out = nc::make_unique<RCSummary>();
  for (const auto& aggregate_and_routes : rc.routes()) {
    const AggregateId& id = aggregate_and_routes.first;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(rc.demands(), id);
    double aggregate_flow_count = demand_and_flow_count.second;

    nc::net::Delay sp_delay = sp.GetDistance(id.src(), id.dst());
    double sp_delay_ms = ToMs(sp_delay);
    sp_delay_ms = std::max(1.0, sp_delay_ms);

    for (const RouteAndFraction& route_and_fraction :
         aggregate_and_routes.second) {
      const nc::net::Walk* path = route_and_fraction.first;
      double fraction = route_and_fraction.second;

      double path_delay_ms = ToMs(path->delay());
      path_delay_ms = std::max(1.0, path_delay_ms);

      double abs_stretch = path_delay_ms - sp_delay_ms;
      double rel_stretch = path_delay_ms / sp_delay_ms;
      CHECK(abs_stretch >= 0);
      CHECK(rel_stretch >= 0);
      CHECK(rel_stretch < 1000);

      size_t flow_count = aggregate_flow_count * fraction;
      flow_count = std::max(1ul, flow_count);
      out->abs_stretches.emplace_back(abs_stretch);
      out->rel_stretches.emplace_back(rel_stretch);
      out->flow_counts.emplace_back(flow_count);
      out->total_sp_delay += sp_delay_ms * flow_count;
      out->total_delay += (sp_delay_ms + abs_stretch) * flow_count;
    }
    out->path_counts.emplace_back(aggregate_and_routes.second.size());
  }

  double locality;
  std::string locality_str =
      nc::FindOrDie(demand_matrix->properties(), "locality");
  CHECK(nc::safe_strtod(locality_str, &locality));

  double load;
  std::string load_str = nc::FindOrDie(demand_matrix->properties(), "load");
  CHECK(nc::safe_strtod(load_str, &load));

  double seed;
  std::string seed_str = nc::FindOrDie(demand_matrix->properties(), "seed");
  CHECK(nc::safe_strtod(seed_str, &seed));

  out->demand_matrix = demand_matrix;
  out->key = std::make_tuple(demand_matrix->graph(), seed, load, locality);

  const nc::net::GraphStorage* graph = rc.graph();
  nc::net::Delay median_link_delay = MedianLinkDelay(*graph);

  nc::net::GraphLinkMap<double> utilizations = rc.LinkUtilizations();
  for (const auto& link_and_utilization : utilizations) {
    nc::net::GraphLinkIndex link = link_and_utilization.first;
    double utilization = *(link_and_utilization.second);
    nc::net::Delay link_delay = graph->GetLink(link)->delay();
    if (link_delay < median_link_delay) {
      out->short_link_utilizations.emplace_back(utilization);
    } else {
      out->long_link_utilizations.emplace_back(utilization);
    }
  }

  return out;
}

// Returns the link scale at which total delay will be 'fraction' of the one at
// no headroom.
static double LinkScaleAtDelayFraction(
    double fraction, const nc::lp::DemandMatrix& demand_matrix) {
  const std::map<std::string, std::string>& properties =
      demand_matrix.properties();
  std::string step_str = nc::FindOrDie(properties, "headroom_vs_delay_step");
  std::string delays_at_scale_str =
      nc::FindOrDie(properties, "headroom_vs_delay");

  double abs_stride;
  CHECK(nc::safe_strtod(step_str, &abs_stride));
  const std::vector<double>& delay_at_scale =
      GetCommaSeparated(delays_at_scale_str);

  std::vector<double> xs;
  std::vector<double> ys;
  for (size_t i = 0; i < delay_at_scale.size(); ++i) {
    double multiplier = 1.0 - i * abs_stride;
    xs.emplace_back(multiplier);
    ys.emplace_back(delay_at_scale[i]);
  }

  nc::Empirical2DFunction f(xs, ys, nc::Empirical2DFunction::LINEAR);
  double min_delay = delay_at_scale.front();
  double max_delay = delay_at_scale.back();
  if (max_delay / min_delay <= fraction) {
    return 1.0 - delay_at_scale.size() * abs_stride;
  }

  for (double link_scale = 1.0; link_scale > 0.0; link_scale -= 0.01) {
    double delay = f.Eval(link_scale);
    double delay_fraction = delay / min_delay;
    if (delay_fraction >= fraction) {
      return link_scale;
    }
  }
  LOG(FATAL) << "Should not happen";
  return 0;
}

class MultiRCSummaryPlotPack {
 public:
  MultiRCSummaryPlotPack()
      : path_stretch_rel_plot_(
            {"Relative path stretch", "fraction longer than SP", "CDF"}),
        path_stretch_max_rel_plot_(
            {"Maximum relative path stretch", "fraction longer than SP"}),
        link_utilization_plot_({"Link utilizations", "utilization"}),
        ratios_plot_({"Total delay", "fraction total delay is longer than SP"}),
        aggregate_path_count_plot_(
            {"Number of paths per aggregate", "path count"}),
        link_scales_plot_(
            {"Link scales at X% increase in total delay", "link scale"}) {}

  void PlotStretchDistribution(const std::string& optimizer,
                               const std::vector<const RCSummary*>& rcs) {
    nc::DiscreteDistribution<uint64_t> dist;
    std::vector<std::pair<double, const RCSummary*>> maxs_decorated;

    for (const RCSummary* rc : rcs) {
      double max = 0;
      for (size_t i = 0; i < rc->rel_stretches.size(); ++i) {
        double rel_stretch = rc->rel_stretches[i];
        double flow_count = rc->flow_counts[i];
        uint64_t value_discrete =
            static_cast<uint64_t>(kDiscreteMultiplier * rel_stretch);
        dist.Add(value_discrete, flow_count);
        max = std::max(rel_stretch, max);
      }
      maxs_decorated.emplace_back(max, rc);
    }

    std::sort(maxs_decorated.begin(), maxs_decorated.end());
    const auto* max = &maxs_decorated.back();
    const auto* med = &maxs_decorated[maxs_decorated.size() / 2];
    path_stretch_max_rel_interesting_.Add("max", max->first, max->second);
    path_stretch_max_rel_interesting_.Add("med", med->first, med->second);

    std::vector<uint64_t> percentiles = dist.Percentiles(kPercentilesCount);
    CHECK(percentiles.size() == kPercentilesCount + 1);

    std::vector<std::pair<double, double>> path_stretch_rel_data;
    for (size_t i = 0; i < percentiles.size(); ++i) {
      double discrete_value = percentiles[i];
      double p = discrete_value / kDiscreteMultiplier;
      path_stretch_rel_data.emplace_back(
          p, static_cast<double>(i) / kPercentilesCount);
    }
    path_stretch_rel_plot_.AddData(optimizer, path_stretch_rel_data);

    std::vector<double> maxs;
    for (const auto& value_and_rc : maxs_decorated) {
      maxs.emplace_back(value_and_rc.first);
    }
    path_stretch_max_rel_plot_.AddData(optimizer, maxs);
  }

  void PlotPathRatios(const std::string& optimizer,
                      const std::vector<const RCSummary*>& rcs) {
    std::vector<std::pair<double, const RCSummary*>> decorated;
    for (const RCSummary* rc : rcs) {
      double change =
          (rc->total_delay - rc->total_sp_delay) / rc->total_sp_delay;
      decorated.emplace_back(change, rc);
    }

    std::sort(decorated.begin(), decorated.end());
    const RCSummary* max = decorated.back().second;
    const RCSummary* med = decorated[decorated.size() / 2].second;

    std::vector<double> ratios;
    for (const auto& change_and_rc : decorated) {
      ratios.emplace_back(change_and_rc.first);
    }

    ratios_interesting_.Add("max", ratios.back(), max);
    ratios_interesting_.Add("med", ratios[ratios.size() / 2], med);
    ratios_plot_.AddData(optimizer, ratios);
  }

  void PlotPathCounts(const std::string& optimizer,
                      const std::vector<const RCSummary*>& rcs) {
    std::vector<std::pair<double, const RCSummary*>> decorated;
    for (const RCSummary* rc : rcs) {
      for (double path_count : rc->path_counts) {
        decorated.emplace_back(path_count, rc);
      }
    }

    std::sort(decorated.begin(), decorated.end());
    const RCSummary* max = decorated.back().second;
    const RCSummary* med = decorated[decorated.size() / 2].second;

    std::vector<double> path_counts;
    for (const auto& value_and_rc : decorated) {
      path_counts.emplace_back(value_and_rc.first);
    }

    aggregate_path_count_interesting_.Add("max", path_counts.back(), max);
    aggregate_path_count_interesting_.Add(
        "med", path_counts[path_counts.size() / 2], med);
    aggregate_path_count_plot_.AddData(optimizer, path_counts);
  }

  void PlotLinkUtilizations(const std::string& optimizer,
                            const std::vector<const RCSummary*>& rcs) {
    std::vector<std::pair<double, const RCSummary*>> decorated;
    for (const RCSummary* rc : rcs) {
      for (double link_utilization : rc->short_link_utilizations) {
        decorated.emplace_back(link_utilization, rc);
      }
      for (double link_utilization : rc->long_link_utilizations) {
        decorated.emplace_back(link_utilization, rc);
      }
    }

    std::sort(decorated.begin(), decorated.end());
    const RCSummary* max = decorated.back().second;
    const RCSummary* med = decorated[decorated.size() / 2].second;

    std::vector<double> link_utilizations;
    for (const auto& value_and_rc : decorated) {
      link_utilizations.emplace_back(value_and_rc.first);
    }

    link_utilization_interesting_.Add("max", link_utilizations.back(), max);
    link_utilization_interesting_.Add(
        "med", link_utilizations[link_utilizations.size() / 2], med);
    link_utilization_plot_.AddData(optimizer, link_utilizations);
  }

  void PlotLinkScalesAtDelay(
      const std::set<const nc::lp::DemandMatrix*>& all_demands) {
    std::vector<double> p_values = {1.01, 1.02, 1.05, 1.10};
    for (double p : p_values) {
      std::vector<double> link_scales;
      for (const nc::lp::DemandMatrix* demand_matrix : all_demands) {
        double scale = LinkScaleAtDelayFraction(p, *demand_matrix);
        link_scales.emplace_back(scale);
      }

      uint32_t percents = (p - 1) * 100;
      link_scales_plot_.AddData(nc::StrCat(percents, "%"), link_scales);
    }
  }

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
};

// Parses a comma-separated string of doubles.
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

class SingleRCSummaryPlotPack {
 public:
  SingleRCSummaryPlotPack()
      : demand_sizes_plot_(
            {"Demand sizes (all possible demands considered)", "size (Mbps)"}),
        sp_utilizations_plot_(
            {"Shortest path link utilization", "link utilization"}),
        cumulative_demands_plot_({"Cumulative distance vs aggregate size",
                                  "SP distance (ms)",
                                  "Cumulative size (Mbps)"}),
        total_delay_at_link_scale_plot_(
            {"Total delay vs link scale", "Link scale", "Total delay"}) {}

  void PlotCumulativeDistances(const nc::lp::DemandMatrix& demand_matrix,
                               const nc::net::AllPairShortestPath& sp) {
    std::vector<std::pair<double, double>> distances_and_demands;
    for (const auto& element : demand_matrix.elements()) {
      nc::net::GraphNodeIndex src = element.src;
      nc::net::GraphNodeIndex dst = element.dst;
      double demand_Mbps = element.demand.Mbps();
      if (demand_Mbps == 0) {
        continue;
      }

      double distance_ms =
          duration_cast<milliseconds>(sp.GetDistance(src, dst)).count();
      distances_and_demands.emplace_back(distance_ms, demand_Mbps);
    }
    std::sort(distances_and_demands.begin(), distances_and_demands.end());

    double total = 0;
    for (auto& distance_and_demand : distances_and_demands) {
      double& demand = distance_and_demand.second;
      total += demand;
      demand = total;
    }

    cumulative_demands_plot_.AddData("", distances_and_demands);
  }

  void PlotDemandSizes(const nc::lp::DemandMatrix& demand_matrix) {
    const nc::net::GraphStorage* graph = demand_matrix.graph();

    std::vector<std::pair<double, ctr::AggregateId>> decorated;
    for (const auto& element : demand_matrix.elements()) {
      nc::net::GraphNodeIndex src = element.src;
      nc::net::GraphNodeIndex dst = element.dst;
      double demand_Mbps = element.demand.Mbps();
      if (demand_Mbps == 0) {
        continue;
      }

      decorated.emplace_back(demand_Mbps, std::make_pair(src, dst));
    }

    std::sort(decorated.begin(), decorated.end());
    ctr::AggregateId max = decorated.back().second;
    ctr::AggregateId med = decorated[decorated.size() / 2].second;

    std::vector<double> demand_sizes;
    for (const auto& value_and_id : decorated) {
      demand_sizes.emplace_back(value_and_id.first);
    }

    demand_sizes_interesting_.Add({"max",
                                   nc::StrCat(demand_sizes.back(), " Mbps"),
                                   max.ToString(*graph)});
    demand_sizes_interesting_.Add(
        {"med", nc::StrCat(demand_sizes[demand_sizes.size() / 2], " Mbps"),
         med.ToString(*graph)});
    demand_sizes_plot_.AddData("", demand_sizes);
  }

  void PlotSPUtilizations(const nc::lp::DemandMatrix& demand_matrix) {
    const nc::net::GraphStorage* graph = demand_matrix.graph();

    std::vector<std::pair<double, nc::net::GraphLinkIndex>> decorated;
    nc::net::GraphLinkMap<double> utilizations = demand_matrix.SPUtilization();
    for (const auto& link_and_utilization : utilizations) {
      decorated.emplace_back(*link_and_utilization.second,
                             link_and_utilization.first);
    }

    std::sort(decorated.begin(), decorated.end());
    nc::net::GraphLinkIndex max = decorated.back().second;
    nc::net::GraphLinkIndex med = decorated[decorated.size() / 2].second;

    std::vector<double> link_utilizations;
    for (const auto& value_and_id : decorated) {
      link_utilizations.emplace_back(value_and_id.first);
    }

    sp_utilizations_interesting_.Add({"max",
                                      nc::StrCat(link_utilizations.back()),
                                      graph->GetLink(max)->ToStringNoPorts()});
    sp_utilizations_interesting_.Add(
        {"med", nc::StrCat(link_utilizations[link_utilizations.size() / 2]),
         graph->GetLink(med)->ToStringNoPorts()});
    sp_utilizations_plot_.AddData("", link_utilizations);
  }

  void PlotTotalDelayAtLinkScale(const nc::lp::DemandMatrix& demand_matrix) {
    const std::map<std::string, std::string>& properties =
        demand_matrix.properties();
    std::string step_str = nc::FindOrDie(properties, "headroom_vs_delay_step");
    std::string delays_at_scale_str =
        nc::FindOrDie(properties, "headroom_vs_delay");

    double abs_stride;
    CHECK(nc::safe_strtod(step_str, &abs_stride));
    const std::vector<double>& delay_at_scale =
        GetCommaSeparated(delays_at_scale_str);

    std::vector<std::pair<double, double>> total_delay_values;
    for (size_t i = 0; i < delay_at_scale.size(); ++i) {
      double multiplier = i * abs_stride;
      total_delay_values.emplace_back(multiplier, delay_at_scale[i]);
    }

    double first_value = total_delay_values.front().second;
    for (auto& multiplier_and_total_delay : total_delay_values) {
      multiplier_and_total_delay.first = 1 - multiplier_and_total_delay.first;
      multiplier_and_total_delay.second =
          multiplier_and_total_delay.second / first_value;
    }
    total_delay_at_link_scale_plot_.AddData("", total_delay_values);
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

  // Plots how adding headroom to all links affects the total delay experienced
  // by all flows in the topology.
  nc::viz::LinePlot total_delay_at_link_scale_plot_;
};

static double FractionOfLinksWithUtilizationMoreThan(
    const std::vector<double>& utilizations, double utilization_threshold) {
  double count = 0;
  for (double utilization : utilizations) {
    if (utilization >= utilization_threshold) {
      ++count;
    }
  }

  return count / utilizations.size();
}

static nc::viz::CDFPlot LinkUtilizationsPlot(
    const std::vector<RCSummary>& rcs,
    const std::vector<double>& utilization_thresholds) {
  nc::viz::CDFPlot out;
  for (double utilization_threshold : utilization_thresholds) {
    std::vector<double> short_utilizations;
    std::vector<double> long_utilizations;
    for (const RCSummary& rc : rcs) {
      double short_f = FractionOfLinksWithUtilizationMoreThan(
          rc.short_link_utilizations, utilization_threshold);
      double long_f = FractionOfLinksWithUtilizationMoreThan(
          rc.long_link_utilizations, utilization_threshold);
      short_utilizations.emplace_back(short_f);
      long_utilizations.emplace_back(long_f);
    }

    out.AddData(nc::StrCat("short >", utilization_threshold, "%"),
                short_utilizations);
    out.AddData(nc::StrCat("long >", utilization_threshold, "%"),
                long_utilizations);
  }

  return out;
}

// static std::string Indent(const std::string& input) {
//  std::vector<std::string> pieces = nc::Split(input, "\n");
//  std::string out;
//  for (const auto& piece : pieces) {
//    nc::StrAppend(&out, "   ", piece, "\n");
//  }
//
//  return out;
//}

class DataPlotter {
 public:
  void AddData(const std::string& opt, const std::string& topology_file,
               std::unique_ptr<RCSummary> rc_summary) {
    const nc::net::GraphStorage* graph;
    double seed;
    double load;
    double locality;
    std::tie(graph, seed, load, locality) = rc_summary->key;

    auto ll = std::make_pair(load, locality);
    graph_to_name_[graph] = topology_file;
    data_by_load_and_locality_[ll][opt].emplace_back(std::move(rc_summary));
  }

  void PlotRoot(const std::string& root) {
    using namespace std::chrono;
    ctemplate::TemplateDictionary dict("plot");

    for (const auto& load_and_locality_and_data : data_by_load_and_locality_) {
      double load;
      double locality;
      std::tie(load, locality) = load_and_locality_and_data.first;

      std::string subdir = nc::StrCat(root, "/data_load_", load, "_", locality);
      const auto& data_map = load_and_locality_and_data.second;
      PlotDataByLoadAndLocality(load, locality, data_map, subdir);

      ctemplate::TemplateDictionary* subdict =
          dict.AddSectionDictionary("subdirs");
      subdict->SetValue("name", nc::StrCat("data_load_", load, "_", locality));
    }

    for (uint32_t i = 0; i < interesting_states.size(); ++i) {
      const RCSummary* rc = interesting_states[i];
      const TMKey& key = rc->key;

      const nc::net::GraphStorage* graph;
      double seed;
      double load;
      double locality;
      std::tie(graph, seed, load, locality) = key;

      const std::string& name = nc::FindOrDie(graph_to_name_, graph);
      std::string subdir = nc::StrCat(root, "/interesting_", i);
      PlotLLMap(name, graph, seed, ll_map, subdir, i);
    }

    std::string output;
    ctemplate::ExpandTemplate(kTopLevelTemplate, ctemplate::DO_NOT_STRIP, &dict,
                              &output);
    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root, "/index.rst"));
  }

 private:
  using LoadAndLocality = std::pair<double, double>;
  using RCSummaryVector = std::vector<std::unique_ptr<RCSummary>>;

  // Plots a single TM. The first argument is summaries, grouped by optimizer.
  // All summaries should be for the same TM.
  void PlotSingleTM(const std::string& root, const RCSummary* rc,
                    const nc::net::AllPairShortestPath& sp,
                    uint32_t interesting_index) {
    using namespace std::chrono;

    const nc::lp::DemandMatrix& demand_matrix = *(rc->demand_matrix);
    SingleRCSummaryPlotPack plot_pack;
    plot_pack.PlotCumulativeDistances(demand_matrix, sp);
    plot_pack.PlotDemandSizes(demand_matrix);
    plot_pack.PlotSPUtilizations(demand_matrix);
    plot_pack.PlotTotalDelayAtLinkScale(demand_matrix);

    ctemplate::TemplateDictionary dict("plot");
    PlotAndAddToTemplate(root, "demand_sizes", plot_pack.demand_sizes_plot_,
                         &dict);
    PlotAndAddToTemplate(root, "sp_utilization",
                         plot_pack.sp_utilizations_plot_, &dict);
    PlotAndAddToTemplate(root, "total_delay_at_link_scale",
                         plot_pack.total_delay_at_link_scale_plot_, &dict);
    PlotAndAddToTemplate(root, "cumulative_demands",
                         plot_pack.cumulative_demands_plot_, &dict);

    AddInterestingDataToTemplate("demand_sizes", {"type", "value", "demand"},
                                 plot_pack.demand_sizes_interesting_, &dict);

    const nc::net::GraphStorage* graph;
    double seed;
    double load;
    double locality;
    std::tie(graph, seed, load, locality) = rc->key;

    dict.SetValue("topology_name", *(rc->topology_file_name));
    dict.SetValue("tm_seed", nc::StrCat(seed));
    dict.SetValue("interesting_index", nc::StrCat(interesting_index));

    std::string output;
    ctemplate::ExpandTemplate(kLocalityLevelTemplate, ctemplate::DO_NOT_STRIP,
                              &dict, &output);
    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root, "/index.rst"));
  }

  void AddInterestingDataToTemplate(const std::string& base,
                                    const std::vector<std::string>& field_ids,
                                    const InterestingTable& data,
                                    ctemplate::TemplateDictionary* dict) {
    const std::vector<std::vector<std::string>>& rows = data.rows();
    const std::vector<const RCSummary*>& rcs = data.rcs();
    bool has_rcs = !rcs.empty();
    if (has_rcs) {
      CHECK(rows.size() == rcs.size());
    }

    for (size_t row_index = 0; row_index < rows.size(); ++row_index) {
      const std::vector<std::string>& row = rows[row_index];
      CHECK(row.size() == field_ids.size());
      ctemplate::TemplateDictionary* subdict =
          dict->AddSectionDictionary(nc::StrCat(base, "_table_rows"));
      for (size_t i = 0; i < field_ids.size(); ++i) {
        const std::string& field_id = field_ids[i];
        const std::string& field_value = row[i];
        subdict->SetValue(field_id, field_value);
      }

      if (has_rcs) {
        const RCSummary* rc = rcs[row_index];
        subdict->SetValue("interesting_index",
                          std::to_string(interesting_data_.size()));
        interesting_data_.emplace_back(rc);
      }
    }
  }

  void PlotAndAddToTemplate(const std::string& root, const std::string& base,
                            const nc::viz::Plot& plot,
                            ctemplate::TemplateDictionary* dict) {
    std::string plot_location = nc::StrCat(root, "/", base, ".svg");
    std::string tgz_location = nc::StrCat(root, "/", base, ".tgz");
    plot.PlotToSVGFile(plot_location);
    plot.PlotToArchiveFile(tgz_location);

    std::string rel_location = nc::StrCat(base, ".svg");
    std::string rel_tgz_location = nc::StrCat(base, ".tgz");
    dict->SetValue(nc::StrCat(base, "_location"), rel_location);
    dict->SetValue(nc::StrCat(base, "_location_tgz"), rel_tgz_location);
  }

  std::map<std::string, std::vector<const RCSummary*>> DataByOptimizer(
      double load, double locality) {
    const std::vector<const RCSummary*>& all_data =
        nc::FindOrDie(data_, std::make_pair(load, locality));

    std::map<std::string, std::vector<const RCSummary*>> out;
    for (const RCSummary* rc : all_data) {
      out[rc->opt].emplace_back(rc);
    }
    return out;
  }

  void PlotDataByLoadAndLocality(const std::string& root, double load,
                                 double locality) {
    MultiRCSummaryPlotPack plot_pack;
    ctemplate::TemplateDictionary dict("plot");
    std::map<std::string, std::vector<const RCSummary*>> opt_map =
        DataByOptimizer(load, locality);

    std::set<const nc::lp::DemandMatrix*> all_demands;
    for (const auto& optimizer_and_rcs : opt_map) {
      const std::string& opt = optimizer_and_rcs.first;
      const std::vector<const RCSummary*>& rcs = optimizer_and_rcs.second;
      for (const RCSummary* rc : rcs) {
        all_demands.emplace(rc->demand_matrix);
      }

      plot_pack.PlotLinkUtilizations(opt, rcs);
      plot_pack.PlotPathCounts(opt, rcs);
      plot_pack.PlotPathRatios(opt, rcs);
      plot_pack.PlotStretchDistribution(opt, rcs);
    }
    plot_pack.PlotLinkScalesAtDelay(all_demands);

    PlotAndAddToTemplate(root, "path_ratios", plots.ratios, &dict);
    PlotAndAddToTemplate(root, "path_stretch_rel", plots.path_stretch_rel,
                         &dict);
    PlotAndAddToTemplate(root, "max_path_stretch_rel",
                         plots.path_stretch_max_rel, &dict);
    PlotAndAddToTemplate(root, "path_count", plots.aggregate_path_count, &dict);
    PlotAndAddToTemplate(root, "link_utilization", plots.link_utilization,
                         &dict);
    PlotAndAddToTemplate(root, "link_scales", link_scales_plot, &dict);
    PlotAndAddToTemplate(root, "ctr_utilizations", ctr_utilizations_plot,
                         &dict);

    dict.SetValue("tm_load", nc::StrCat(load));
    dict.SetValue("tm_locality", nc::StrCat(locality));
    std::string output;
    ctemplate::ExpandTemplate(kSummaryTemplate, ctemplate::DO_NOT_STRIP, &dict,
                              &output);
    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root, "/index.rst"));
  }

  // All data, grouped by load and locality.
  std::map<LoadAndLocality, std::vector<const RCSummary*>> data_;

  std::vector<const RCSummary*> interesting_data_;
};

static std::string GetFilename(const std::string& tm_file,
                               const std::string opt_string) {
  std::string tm_base = nc::StringReplace(tm_file, ".demands", "", true);
  return nc::StrCat(tm_base, "_", opt_string, ".rc");
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::tie(topologies, matrices) = ctr::GetDemandMatrixInputs(false);
  std::vector<std::string> optimizers = nc::Split(FLAGS_optimizers, ",");

  // Need to know the shortest paths for each graph.
  std::map<const nc::net::GraphStorage*, nc::net::AllPairShortestPath> sp_map;

  std::map<std::string, const ctr::TopologyAndFilename*> topologies_by_name;
  for (const auto& topology : topologies) {
    topologies_by_name[topology.file] = &topology;

    const nc::net::GraphStorage* graph = topology.graph.get();
    sp_map.emplace(
        std::piecewise_construct, std::forward_as_tuple(graph),
        std::forward_as_tuple(nc::net::ExclusionSet(), graph->AdjacencyList(),
                              nullptr, nullptr));
  }

  // A tree of plots.
  DataPlotter data_plotter;
  for (const auto& matrix : matrices) {
    const ctr::TopologyAndFilename* topology_and_filename =
        nc::FindOrDieNoPrint(topologies_by_name, matrix.topology_file);

    ctr::PathProvider path_provider(topology_and_filename->graph.get());
    const nc::lp::DemandMatrix* demand_matrix = matrix.demand_matrix.get();

    std::vector<std::unique_ptr<RoutingConfiguration>> rcs;
    for (const std::string& opt : optimizers) {
      std::string rc_filename = GetFilename(matrix.file, opt);
      if (!nc::File::Exists(rc_filename)) {
        break;
      }

      std::string rc_serialized = nc::File::ReadFileToStringOrDie(rc_filename);
      std::unique_ptr<TrafficMatrix> tm =
          TrafficMatrix::DistributeFromDemandMatrix(*demand_matrix);
      LOG(INFO) << "Will parse " << matrix.file << " at "
                << matrix.topology_file << " opt " << opt;
      auto rc = RoutingConfiguration::LoadFromSerializedText(
          *tm, topology_and_filename->node_order, rc_serialized,
          &path_provider);
      rcs.emplace_back(std::move(rc));
    }

    if (rcs.size() != optimizers.size()) {
      LOG(INFO) << "Will skip " << matrix.file << " at "
                << matrix.topology_file;
      continue;
    }

    for (size_t i = 0; i < optimizers.size(); ++i) {
      const std::string& opt = optimizers[i];
      const RoutingConfiguration& rc = *(rcs[i]);

      const nc::net::AllPairShortestPath& sp =
          nc::FindOrDieNoPrint(sp_map, rc.graph());
      auto summary = ParseRC(rc, sp, demand_matrix);
      data_plotter.AddData(opt, matrix.topology_file, std::move(summary));
    }
  }

  data_plotter.PlotRoot("plot_tree");
}
