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

static RCSummary ParseRC(const RoutingConfiguration& rc,
                         const nc::net::AllPairShortestPath& sp,
                         const nc::lp::DemandMatrix* demand_matrix) {
  RCSummary out;
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
      double rel_stretch = (path_delay_ms - sp_delay_ms) / sp_delay_ms;
      CHECK(abs_stretch >= 0);
      CHECK(rel_stretch >= 0);
      CHECK(rel_stretch < 1000);

      size_t flow_count = aggregate_flow_count * fraction;
      flow_count = std::max(1ul, flow_count);
      out.abs_stretches.emplace_back(abs_stretch);
      out.rel_stretches.emplace_back(rel_stretch);
      out.flow_counts.emplace_back(flow_count);
      out.total_sp_delay += sp_delay_ms * flow_count;
      out.total_delay += (sp_delay_ms + abs_stretch) * flow_count;
    }
    out.path_counts.emplace_back(aggregate_and_routes.second.size());
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

  out.demand_matrix = demand_matrix;
  out.key = std::make_tuple(demand_matrix->graph(), seed, load, locality);

  const nc::net::GraphStorage* graph = rc.graph();
  nc::net::Delay median_link_delay = MedianLinkDelay(*graph);

  nc::net::GraphLinkMap<double> utilizations = rc.LinkUtilizations();
  for (const auto& link_and_utilization : utilizations) {
    nc::net::GraphLinkIndex link = link_and_utilization.first;
    double utilization = *(link_and_utilization.second);
    nc::net::Delay link_delay = graph->GetLink(link)->delay();
    if (link_delay < median_link_delay) {
      out.short_link_utilizations.emplace_back(utilization);
    } else {
      out.long_link_utilizations.emplace_back(utilization);
    }
  }

  return out;
}

static std::tuple<std::vector<double>, const RCSummary*, const RCSummary*,
                  double, double>
GetPathRatios(const std::vector<RCSummary>& rcs) {
  std::vector<std::pair<double, const RCSummary*>> decorated;
  for (const RCSummary& rc : rcs) {
    double change = (rc.total_delay - rc.total_sp_delay) / rc.total_sp_delay;
    decorated.emplace_back(change, &rc);
  }

  std::sort(decorated.begin(), decorated.end());
  const RCSummary* max = decorated.back().second;
  const RCSummary* med = decorated[decorated.size() / 2].second;

  std::vector<double> out;
  for (const auto& change_and_state : decorated) {
    out.emplace_back(change_and_state.first);
  }

  return std::make_tuple(out, med, max, out[out.size() / 2], out.back());
}

static std::pair<std::vector<double>, std::vector<double>>
GetStretchDistribution(const std::vector<RCSummary>& rcs) {
  nc::DiscreteDistribution<uint64_t> dist;
  std::vector<double> maxs;

  for (const RCSummary& rc : rcs) {
    double max = 0;
    for (size_t i = 0; i < rc.rel_stretches.size(); ++i) {
      double rel_stretch = rc.rel_stretches[i];
      double flow_count = rc.flow_counts[i];
      uint64_t value_discrete =
          static_cast<uint64_t>(kDiscreteMultiplier * rel_stretch);
      dist.Add(value_discrete, flow_count);
      max = std::max(rel_stretch, max);
    }
    maxs.emplace_back(max);
  }

  std::vector<uint64_t> percentiles = dist.Percentiles(kPercentilesCount);
  CHECK(percentiles.size() == kPercentilesCount + 1);

  std::vector<double> out;
  out.reserve(percentiles.size());
  for (uint64_t discrete_value : percentiles) {
    double v = static_cast<double>(discrete_value) / kDiscreteMultiplier;
    out.emplace_back(v);
  }

  return {out, maxs};
}

struct PlotPack {
  PlotPack()
      : path_stretch_rel(
            {"Relative path stretch", "fraction longer than SP", "CDF"}),
        path_stretch_max_rel(
            {"Maximum relative path stretch", "fraction longer than SP"}),
        link_utilization({"Link utilizations", "utilization"}),
        ratios({"Total delay", "fraction total delay is longer than SP"}),
        aggregate_path_count({"Number of paths per aggregate", "path count"}) {}

  nc::viz::LinePlot path_stretch_rel;
  nc::viz::CDFPlot path_stretch_max_rel;
  nc::viz::CDFPlot link_utilization;
  nc::viz::CDFPlot ratios;
  nc::viz::CDFPlot aggregate_path_count;
};

struct HandleSingleOptimizerResult {
  const RCSummary* med;
  const RCSummary* max;
  double med_value;
  double max_value;
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

static HandleSingleOptimizerResult HandleSingleOptimizer(
    const std::string& optimizer, const std::vector<RCSummary>& rcs,
    PlotPack* plots) {
  std::vector<double> ratios;
  const RCSummary* med;
  const RCSummary* max;
  double med_value;
  double max_value;
  std::tie(ratios, med, max, med_value, max_value) = GetPathRatios(rcs);
  plots->ratios.AddData(optimizer, ratios);
  CHECK(med != nullptr);
  CHECK(max != nullptr);

  // Will plot the percentiles to get a CDF.
  std::vector<double> percentiles;
  std::vector<double> maxs;
  std::tie(percentiles, maxs) = GetStretchDistribution(rcs);
  std::vector<std::pair<double, double>> path_stretch_rel_data;
  for (size_t i = 0; i < percentiles.size(); ++i) {
    double p = percentiles[i];
    path_stretch_rel_data.emplace_back(p, i);
  }
  plots->path_stretch_rel.AddData(optimizer, path_stretch_rel_data);
  plots->path_stretch_max_rel.AddData(optimizer, maxs);

  std::vector<double> per_aggregate_path_counts;
  for (const RCSummary& rc : rcs) {
    per_aggregate_path_counts.insert(per_aggregate_path_counts.end(),
                                     rc.path_counts.begin(),
                                     rc.path_counts.end());
  }
  plots->aggregate_path_count.AddData(optimizer, per_aggregate_path_counts);

  std::vector<double> link_utilizations;
  for (const RCSummary& rc : rcs) {
    link_utilizations.insert(link_utilizations.end(),
                             rc.short_link_utilizations.begin(),
                             rc.short_link_utilizations.end());
    link_utilizations.insert(link_utilizations.end(),
                             rc.long_link_utilizations.begin(),
                             rc.long_link_utilizations.end());
  }
  plots->link_utilization.AddData(optimizer, link_utilizations);

  return {med, max, med_value, max_value};
}

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
               const RCSummary& rc_summary) {
    const nc::net::GraphStorage* graph;
    double seed;
    double load;
    double locality;
    std::tie(graph, seed, load, locality) = rc_summary.key;

    auto ll = std::make_pair(load, locality);
    graph_to_name_[graph] = topology_file;
    top_level_map_[graph][seed][ll] = rc_summary.demand_matrix;
    data_by_load_and_locality_[ll][opt].emplace_back(rc_summary);
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
      const TMKey& key = interesting_states[i];

      const nc::net::GraphStorage* graph;
      double seed;
      double load;
      double locality;
      std::tie(graph, seed, load, locality) = key;

      const std::string& name = nc::FindOrDie(graph_to_name_, graph);
      const LLMap& ll_map =
          nc::FindOrDie(nc::FindOrDie(top_level_map_, graph), seed);
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

  void PlotSingleTM(const RCSummary& summary, uint32_t interesting_index) {
    using namespace std::chrono;

    nc::viz::CDFPlot demand_sizes_plot(
        {"Demand sizes (all possible demands considered)", "size (Mbps)"});
    nc::viz::CDFPlot sp_utilizations_plot(
        {"Shortest path link utilization", "link utilization"});
    nc::viz::LinePlot cumulative_demands_plot(
        {"Cumulative distance vs aggregate size", "SP distance (ms)",
         "Cumulative size (Mbps)"});
    nc::viz::LinePlot total_delay_at_link_scale_plot(
        {"Total delay vs link scale", "Link scale", "Total delay"});

    ctemplate::TemplateDictionary dict("plot");
    for (const auto& load_locality_and_tm : ll_map) {
      double load;
      double locality;
      std::tie(load, locality) = load_locality_and_tm.first;
      const nc::lp::DemandMatrix* demand_matrix = load_locality_and_tm.second;

      std::vector<double> demand_sizes_Mbps;
      std::vector<double> sp_utilizations;
      std::vector<std::pair<double, double>> cumulative_distances;

      nc::net::AllPairShortestPath sp({}, graph->AdjacencyList(), nullptr,
                                      nullptr);
      for (const auto& element : demand_matrix->elements()) {
        double demand_Mbps = element.demand.Mbps();
        if (demand_Mbps == 0) {
          continue;
        }

        double distance_ms = duration_cast<milliseconds>(
                                 sp.GetDistance(element.src, element.dst))
                                 .count();
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

      nc::net::GraphLinkMap<double> utilizations =
          demand_matrix->SPUtilization();
      for (const auto& link_and_utilization : utilizations) {
        sp_utilizations.emplace_back(*link_and_utilization.second);
      }

      size_t node_count = graph->NodeCount();
      double possible_count = node_count * (node_count - 1);
      CHECK(demand_sizes_Mbps.size() <= possible_count);
      std::string id = nc::StrCat("load ", load, ", locality ", locality);

      demand_sizes_plot.AddData(id, demand_sizes_Mbps);
      sp_utilizations_plot.AddData(id, sp_utilizations);
      cumulative_demands_plot.AddData(id, cumulative_distances);

      const std::map<std::string, std::string>& properties =
          demand_matrix->properties();
      std::string step_str =
          nc::FindOrDie(properties, "headroom_vs_delay_step");
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
      total_delay_at_link_scale_plot.AddData(id, total_delay_values);

      ctemplate::TemplateDictionary* subdict =
          dict.AddSectionDictionary("summary_table_rows");

      size_t element_count = demand_matrix->elements().size();
      double element_count_fraction = element_count / possible_count;
      subdict->SetValue("load_value", nc::ToStringMaxDecimals(load, 2));
      subdict->SetValue("locality_value", nc::ToStringMaxDecimals(locality, 2));
      subdict->SetValue("demand_count_value", std::to_string(element_count));
      subdict->SetValue("demand_fraction_value",
                        nc::ToStringMaxDecimals(element_count_fraction, 3));
    }

    PlotAndAddToTemplate(root, "demand_sizes", demand_sizes_plot, &dict);
    PlotAndAddToTemplate(root, "sp_utilization", sp_utilizations_plot, &dict);
    PlotAndAddToTemplate(root, "total_delay_at_link_scale",
                         total_delay_at_link_scale_plot, &dict);
    PlotAndAddToTemplate(root, "cumulative_demands", cumulative_demands_plot,
                         &dict);

    dict.SetValue("topology_name", topology_name);
    dict.SetValue("tm_seed", nc::StrCat(seed));
    dict.SetValue("interesting_index", nc::StrCat(interesting_index));

    std::string output;
    ctemplate::ExpandTemplate(kLocalityLevelTemplate, ctemplate::DO_NOT_STRIP,
                              &dict, &output);
    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root, "/index.rst"));
  }

  //  void PlotSeedMap(const std::string& name, const nc::net::GraphStorage*
  //  graph,
  //                   const SeedMap& seed_map, const std::string& root) {
  //    using namespace std::chrono;
  //    nc::viz::CDFPlot node_distances_plot(
  //        {"Length of each of all N*(N - 1) pair\\'s shortest path",
  //         "SP length (ms)"});
  //
  //    nc::net::AllPairShortestPath sp({}, graph->AdjacencyList(), nullptr,
  //                                    nullptr);
  //    std::vector<double> node_distances_ms;
  //    for (nc::net::GraphNodeIndex src : graph->AllNodes()) {
  //      for (nc::net::GraphNodeIndex dst : graph->AllNodes()) {
  //        if (src == dst) {
  //          continue;
  //        }
  //
  //        double distance_ms =
  //            duration_cast<milliseconds>(sp.GetDistance(src, dst)).count();
  //        node_distances_ms.emplace_back(distance_ms);
  //      }
  //    }
  //    node_distances_plot.AddData("", node_distances_ms);
  //
  //    ctemplate::TemplateDictionary dict("plot");
  //    dict.SetValue("topology_name", name);
  //
  //    for (const auto& seed_and_ll : seed_map) {
  //      double seed = seed_and_ll.first;
  //      const LLMap& ll_map = seed_and_ll.second;
  //      std::string subdir = nc::StrCat(root, "/tm_seed_", seed);
  //      PlotLLMap(name, graph, seed, ll_map, subdir);
  //
  //      ctemplate::TemplateDictionary* subdict =
  //          dict.AddSectionDictionary("traffic_matrices");
  //      subdict->SetValue("tm_seed", nc::StrCat(seed));
  //    }
  //
  //    nc::net::GraphStats stats = graph->Stats();
  //    dict.SetValue("topology_summary", Indent(stats.ToString()));
  //
  //    PlotAndAddToTemplate(root, "node_distances", node_distances_plot,
  //    &dict);
  //
  //    std::string output;
  //    ctemplate::ExpandTemplate(kSingleTopTemplate, ctemplate::DO_NOT_STRIP,
  //                              &dict, &output);
  //    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root,
  //    "/index.rst"));
  //  }

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

  void PlotDataByLoadAndLocality(
      double load, double locality,
      const std::map<std::string, std::vector<RCSummary>>& opt_map,
      const std::string& root) {
    PlotPack plots;
    ctemplate::TemplateDictionary dict("plot");
    std::set<const nc::lp::DemandMatrix*> all_demands;

    for (const auto& optimizer_and_rcs : opt_map) {
      const std::string& opt = optimizer_and_rcs.first;
      const std::vector<RCSummary>& rcs = optimizer_and_rcs.second;
      for (const RCSummary& rc : rcs) {
        all_demands.emplace(rc.demand_matrix);
      }

      HandleSingleOptimizerResult result =
          HandleSingleOptimizer(opt, rcs, &plots);
      ctemplate::TemplateDictionary* subdict =
          dict.AddSectionDictionary("max_ratio_table_rows");
      subdict->SetValue("optimizer", opt);
      subdict->SetValue("value", nc::StrCat(result.max_value));
      subdict->SetValue("interesting_index",
                        std::to_string(interesting_states.size()));
      interesting_states.emplace_back(result.max->key);

      subdict = dict.AddSectionDictionary("med_ratio_table_rows");
      subdict->SetValue("optimizer", opt);
      subdict->SetValue("value", nc::StrCat(result.med_value));
      subdict->SetValue("interesting_index",
                        std::to_string(interesting_states.size()));
      interesting_states.emplace_back(result.med->key);
    }

    nc::viz::CDFPlot link_scales_plot(
        {"Link scales at X% increase in total delay", "link scale"});
    std::vector<double> p_values = {1.01, 1.02, 1.05, 1.10};
    for (double p : p_values) {
      std::vector<double> link_scales;
      for (const nc::lp::DemandMatrix* demand_matrix : all_demands) {
        double scale = LinkScaleAtDelayFraction(p, *demand_matrix);
        link_scales.emplace_back(scale);
      }

      uint32_t percents = (p - 1) * 100;
      link_scales_plot.AddData(nc::StrCat(percents, "%"), link_scales);
    }

    // Will only plot link utilizations for CTR.
    const std::vector<RCSummary>& ctr_rcs = nc::FindOrDie(opt_map, "CTR");
    nc::viz::CDFPlot ctr_utilizations_plot =
        LinkUtilizationsPlot(ctr_rcs, {0.9, 0.6});

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

  // Stores graphs.
  std::map<const nc::net::GraphStorage*, std::string> graph_to_name_;

  // All data, grouped by load and locality.
  std::map<LoadAndLocality, std::map<std::string, std::vector<RCSummary>>>
      data_by_load_and_locality_;

  // Stores all data.
  TopologyMap top_level_map_;

  // Interesting states to display.
  std::vector<TMKey> interesting_states;
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
      RCSummary summary = ParseRC(rc, sp, demand_matrix);
      data_plotter.AddData(opt, matrix.topology_file, summary);
    }
  }

  data_plotter.PlotRoot("plot_tree");
}
