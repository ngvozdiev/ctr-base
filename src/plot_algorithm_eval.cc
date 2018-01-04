#include <gflags/gflags.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode/common.h"
#include "ncode/file.h"
#include "ncode/logging.h"
#include "ncode/lp/demand_matrix.h"
#include "ncode/map_util.h"
#include "ncode/net/algorithm.h"
#include "ncode/net/net_common.h"
#include "ncode/net/net_gen.h"
#include "ncode/strutil.h"
#include "ncode/viz/grapher.h"
#include "ncode/viz/ctemplate/template.h"
#include "ncode/viz/ctemplate/template_dictionary.h"
#include "ncode/viz/ctemplate/template_enums.h"

DEFINE_string(metrics_dir, "", "The metrics directory.");
DEFINE_string(optimizers, "MinMaxLD,CTRNFC,B4,CTR,MinMaxK10",
              "Optimizers to plot.");
DEFINE_string(output_prefix, "",
              "Prefix that will be added to all output directories.");

static constexpr char kTopLevelTemplate[] = "top_level_template.rst";
static constexpr char kSingleTopTemplate[] = "single_top_template.rst";
static constexpr char kSummaryTemplate[] = "summary_template.rst";
static constexpr char kLocalityLevelTemplate[] = "tm_stat_template.rst";
static constexpr char kSingleTmTemplate[] = "single_tm_template.rst";

static constexpr size_t kDiscreteMultiplier = 10000;
static constexpr size_t kPercentilesCount = 10000;

using namespace nc::metrics::parser;
using namespace ctr::alg_eval;
using namespace ctr;
using namespace std::chrono;

using TopologyAndTM = std::pair<std::string, std::string>;

static double ToMs(nc::net::Delay delay) {
  return duration_cast<microseconds>(delay).count() / 1000.0;
}

static std::tuple<std::vector<double>, const RoutingConfiguration*,
                  const RoutingConfiguration*, double, double>
GetPathRatios(const std::vector<const RoutingConfiguration*>& routing_configs) {
  std::vector<std::pair<double, const RoutingConfiguration*>> decorated;
  for (const RoutingConfiguration* rc : routing_configs) {
    double total_sp_delay = 0;
    double total_delay = 0;
    for (const auto& aggregate_and_routes : rc->routes()) {
      const AggregateId& id = aggregate_and_routes.first;
      const DemandAndFlowCount& demand_and_flow_count =
          nc::FindOrDie(rc->demands(), id);
      double aggregate_flow_count = demand_and_flow_count.second;

      nc::net::Delay sp_delay = id.GetSPDelay(*rc->graph());
      double sp_delay_ms = ToMs(sp_delay);

      for (const RouteAndFraction& route_and_fraction :
           aggregate_and_routes.second) {
        const nc::net::Walk* path = route_and_fraction.first;
        double fraction = route_and_fraction.second;

        double path_delay_ms = ToMs(path->delay());
        double abs_stretch = path_delay_ms - sp_delay_ms;
        size_t flow_count = aggregate_flow_count * fraction;
        flow_count = std::max(1ul, flow_count);

        total_sp_delay += sp_delay_ms * flow_count;
        total_delay += (sp_delay_ms + abs_stretch) * flow_count;
      }
    }

    double change = (total_delay - total_sp_delay) / total_sp_delay;
    decorated.emplace_back(change, rc);
  }

  std::sort(decorated.begin(), decorated.end());
  const RoutingConfiguration* max = decorated.back().second;
  const RoutingConfiguration* med = decorated[decorated.size() / 2].second;

  std::vector<double> out;
  for (const auto& change_and_state : decorated) {
    out.emplace_back(change_and_state.first);
  }

  return std::make_tuple(out, med, max, out[out.size() / 2], out.back());
}

static std::pair<std::vector<double>, std::vector<double>>
GetStretchDistribution(const std::vector<const RoutingConfiguration*>& rcs) {
  nc::DiscreteDistribution<uint64_t> dist;
  std::vector<double> maxs;

  for (const RoutingConfiguration* rc : rcs) {
    double max = 0;
    for (const auto& aggregate_and_routes : rc->routes()) {
      const AggregateId& id = aggregate_and_routes.first;
      const DemandAndFlowCount& demand_and_flow_count =
          nc::FindOrDie(rc->demands(), id);
      double aggregate_flow_count = demand_and_flow_count.second;

      nc::net::Delay sp_delay = id.GetSPDelay(*rc->graph());
      double sp_delay_ms = ToMs(sp_delay);

      for (const RouteAndFraction& route_and_fraction :
           aggregate_and_routes.second) {
        const nc::net::Walk* path = route_and_fraction.first;
        double fraction = route_and_fraction.second;

        double path_delay_ms = ToMs(path->delay());
        size_t flow_count = aggregate_flow_count * fraction;

        CHECK(path_delay_ms >= sp_delay_ms);
        double abs_stretch = path_delay_ms - sp_delay_ms;
        double rel_stretch = abs_stretch / sp_delay_ms;
        max = std::max(rel_stretch, max);

        // Have to discretize the value.
        uint64_t value_discrete =
            static_cast<uint64_t>(kDiscreteMultiplier * abs_stretch);
        dist.Add(value_discrete, flow_count);
      }
    }
    maxs.emplace_back(max);
  }

  std::vector<uint64_t> percentiles = dist.Percentiles(kPercentilesCount);
  CHECK(percentiles.size() == kPercentilesCount + 1);

  std::vector<double> out;
  out.reserve(percentiles.size());
  for (uint64_t discrete_value : percentiles) {
    out.emplace_back(static_cast<double>(discrete_value) / kDiscreteMultiplier);
  }

  return {out, maxs};
}

// static std::vector<double> GetTotalDelays(const TMStateMap& tm_state_map) {
//  std::vector<double> out;
//  for (const auto& key_and_data : tm_state_map) {
//    const OptimizerTMState& tm_state = *(key_and_data.second);
//    std::vector<AggregateTMState> aggregates = tm_state.GetAggregates();
//    double total_delay = 0;
//    for (const auto& aggregate : aggregates) {
//      double sp_delay_ms = aggregate.sp_delay_ms;
//      for (const auto& delay_and_count : aggregate.paths) {
//        double path_delay_ms = delay_and_count.first;
//        double abs_stretch = path_delay_ms - sp_delay_ms;
//        uint32_t flow_count = delay_and_count.second;
//        total_delay += (sp_delay_ms + abs_stretch) * flow_count;
//      }
//    }
//
//    out.emplace_back(total_delay);
//  }
//
//  return out;
//}

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
  const RoutingConfiguration* med;
  const RoutingConfiguration* max;
  double med_value;
  double max_value;
};

static HandleSingleOptimizerResult HandleSingleOptimizer(
    const std::string& optimizer,
    const std::vector<const RoutingConfiguration*>& rcs, PlotPack* plots) {
  std::vector<double> ratios;
  const RoutingConfiguration* med;
  const RoutingConfiguration* max;
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
  for (const RoutingConfiguration* rc : rcs) {
    std::vector<size_t> num_paths = rc->NumberOfPathsInAggregate();
    for (size_t num : num_paths) {
      per_aggregate_path_counts.emplace_back(num);
    }
  }
  plots->aggregate_path_count.AddData(optimizer, per_aggregate_path_counts);

  std::vector<double> link_utilizations;
  for (const RoutingConfiguration* rc : rcs) {
    const nc::net::GraphLinkMap<double>& utilizations = rc->LinkUtilizations();
    for (const auto& link_and_utilization : utilizations) {
      link_utilizations.emplace_back(*link_and_utilization.second);
    }
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

// static void PlotTotalDelay(const std::map<std::string, TMStateMap>& data,
//                           const std::vector<std::string>& optimizers) {
//  std::vector<nc::viz::DataSeries2D> to_plot;
//  for (const std::string& opt : optimizers) {
//    const TMStateMap& state_map = nc::FindOrDie(data, opt);
//
//    std::vector<std::pair<double, double>> xy;
//    std::vector<double> total_delays = GetTotalDelays(state_map);
//    for (size_t i = 0; i < total_delays.size(); ++i) {
//      xy.emplace_back(i, total_delays[i]);
//    }
//
//    to_plot.push_back({opt, xy});
//  }
//
//  nc::viz::LinePlot plot;
//  plot.AddData(to_plot);
//  plot.PlotToDir("total_delays");
//}

// Returns a list of topology/tm pairs ordered by number of aggregates.
// static std::vector<std::pair<size_t, TopologyAndTM>>
// TopologiesAndTMOrderedBySize(const TMStateMap& data) {
//  std::vector<std::pair<size_t, TopologyAndTM>> out;
//  for (const auto& topology_and_tm_and_state : data) {
//    const TopologyAndTM& topology_and_tm = topology_and_tm_and_state.first;
//    const OptimizerTMState& tm_state = *(topology_and_tm_and_state.second);
//
//    size_t aggregate_count = tm_state.aggregate_path_count.size();
//    out.emplace_back(aggregate_count, topology_and_tm);
//  }
//
//  std::sort(out.begin(), out.end());
//  return out;
//}
//
// static void PlotCTRRuntime(const std::map<std::string, TMStateMap>& data) {
//  DataMap runtime_data = SimpleParseNumericDataNoTimestamps(
//      FLAGS_metrics_dir, kRuntimeMetric, ".*");
//  DataMap runtime_data_cached = SimpleParseNumericDataNoTimestamps(
//      FLAGS_metrics_dir, kRuntimeCachedMetric, ".*");
//
//  const TMStateMap& state_map = nc::FindOrDie(data, "CTR");
//  std::vector<std::pair<size_t, TopologyAndTM>> runs_ordered =
//      TopologiesAndTMOrderedBySize(state_map);
//
//  std::vector<double> y;
//  std::vector<double> y_cached;
//  for (const auto& size_and_topology : runs_ordered) {
//    size_t size = size_and_topology.first;
//    const TopologyAndTM& topology_and_tm = size_and_topology.second;
//
//    std::string key =
//        nc::StrCat(topology_and_tm.first, ":", topology_and_tm.second);
//    const DataVector& runtime_data_vector =
//        nc::FindOrDieNoPrint(runtime_data, std::make_pair(kRuntimeMetric,
//        key));
//    const DataVector& runtime_data_cached_vector = nc::FindOrDieNoPrint(
//        runtime_data_cached, std::make_pair(kRuntimeCachedMetric, key));
//
//    CHECK(runtime_data_vector.size() == 1);
//    CHECK(runtime_data_cached_vector.size() == 1);
//
//    y.emplace_back(runtime_data_vector.front());
//    y_cached.emplace_back(runtime_data_cached_vector.front());
//  }
//
//  nc::viz::CDFPlot plot;
//  plot.AddData("Regular", y_cached);
//  plot.AddData("Cold run", y);
//  plot.PlotToDir("runtime");
//}

static std::string Indent(const std::string& input) {
  std::vector<std::string> pieces = nc::Split(input, "\n");
  std::string out;
  for (const auto& piece : pieces) {
    nc::StrAppend(&out, "   ", piece, "\n");
  }

  return out;
}

class DataPlotter {
 public:
  void AddData(const std::string& opt, const TopologyAndTM& topology_and_tm,
               const RoutingConfiguration* rc) {
    const std::string& topology_file = topology_and_tm.first;
    const std::string& tm_file = topology_and_tm.second;
    std::unique_ptr<GraphAndNodeOrder>& graph_and_node_order =
        graphs_[topology_file];

    // In case the same graph was already used by a different TM.
    if (!graph_and_node_order) {
      std::vector<std::string> new_node_order;
      nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
          nc::File::ReadFileToStringOrDie(topology_file), &new_node_order);
      builder.RemoveMultipleLinks();
      graph_and_node_order =
          nc::make_unique<GraphAndNodeOrder>(builder, new_node_order);

      std::string name = nc::File::ExtractFileName(topology_file);
      name = nc::StringReplace(name, ".graph", "", true);

      graph_to_name_[&graph_and_node_order->graph] = name;
    }

    // Same as above, but for TMs.
    const nc::net::GraphStorage* graph = &graph_and_node_order->graph;
    std::unique_ptr<nc::lp::DemandMatrix>& demand_matrix =
        demands_by_topology_[topology_and_tm];
    if (!demand_matrix) {
      const std::vector<std::string>& node_order =
          graph_and_node_order->node_order;
      demand_matrix = nc::lp::DemandMatrix::LoadRepetitaFileOrDie(
          tm_file, node_order, graph);
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

    TMKey key = std::make_tuple(graph, seed, load, locality);
    demands_by_key_[key] = demand_matrix.get();
    key_by_state_[rc] = key;

    auto ll = std::make_pair(load, locality);
    top_level_map_[graph][seed][ll] = demand_matrix.get();
    data_by_load_and_localiy_[ll][opt].emplace_back(rc);
  }

  void PlotRoot(const std::string& root) {
    using namespace std::chrono;
    ctemplate::TemplateDictionary dict("plot");

    for (const auto& load_and_locality_and_data : data_by_load_and_localiy_) {
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

    //    for (const auto& graph_and_seed_map : top_level_map_) {
    //      const nc::net::GraphStorage* graph = graph_and_seed_map.first;
    //      const SeedMap& seed_map = graph_and_seed_map.second;
    //      const std::string& name = nc::FindOrDie(graph_to_name_, graph);
    //
    //      std::string subdir = nc::StrCat(root, "/", name);
    //      PlotSeedMap(name, graph, seed_map, subdir);
    //
    //      nc::net::GraphStats stats = graph->Stats();
    //      double diameter_ms =
    //          duration_cast<milliseconds>(stats.sp_delay_percentiles.back())
    //              .count();
    //      double min_linkspeed =
    //      stats.link_capacity_percentiles.front().Mbps();
    //      double max_linkspeed =
    //      stats.link_capacity_percentiles.back().Mbps();
    //
    //      ctemplate::TemplateDictionary* subdict =
    //          dict.AddSectionDictionary("summary_table_rows");
    //      subdict->SetValue("name", name);
    //      subdict->SetValue("node_count", std::to_string(stats.nodes_count));
    //      subdict->SetValue("link_count", std::to_string(stats.links_count));
    //      subdict->SetValue("diameter", nc::StrCat(diameter_ms, " ms"));
    //      subdict->SetValue("linkspeed", nc::StrCat(min_linkspeed, " Mbps / ",
    //                                                max_linkspeed, " Mbps"));
    //
    //      subdict = dict.AddSectionDictionary("subdirs");
    //      subdict->SetValue("name", name);
    //    }

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
      PlotLLMap(name, graph, seed, ll_map, subdir);
    }

    std::string output;
    ctemplate::ExpandTemplate(kTopLevelTemplate, ctemplate::DO_NOT_STRIP, &dict,
                              &output);
    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root, "/main.rst"));
  }

 private:
  // Graph, seed and load and locality.
  using TMKey =
      std::tuple<const nc::net::GraphStorage*, double, double, double>;
  using LoadAndLocality = std::pair<double, double>;

  struct GraphAndNodeOrder {
    GraphAndNodeOrder(const nc::net::GraphBuilder& builder,
                      const std::vector<std::string>& node_order)
        : graph(builder), node_order(node_order) {}

    nc::net::GraphStorage graph;
    std::vector<std::string> node_order;
  };

  // Maps a combination of load and locality to a traffic matrix.
  using LLMap = std::map<LoadAndLocality, const nc::lp::DemandMatrix*>;

  // Groups TMs with the same seed value.
  using SeedMap = std::map<double, LLMap>;

  // Maps a topology to all traffic matrices on the same topology.
  using TopologyMap = std::map<const nc::net::GraphStorage*, SeedMap>;

  void PlotLLMap(const std::string& topology_name,
                 const nc::net::GraphStorage* graph, double seed,
                 const LLMap& ll_map, const std::string& root) {
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

    std::string output;
    ctemplate::ExpandTemplate(kLocalityLevelTemplate, ctemplate::DO_NOT_STRIP,
                              &dict, &output);
    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root, "/main.rst"));
  }

  void PlotSeedMap(const std::string& name, const nc::net::GraphStorage* graph,
                   const SeedMap& seed_map, const std::string& root) {
    using namespace std::chrono;
    nc::viz::CDFPlot node_distances_plot(
        {"Length of each of all N*(N - 1) pair\\'s shortest path",
         "SP length (ms)"});

    nc::net::AllPairShortestPath sp({}, graph->AdjacencyList(), nullptr,
                                    nullptr);
    std::vector<double> node_distances_ms;
    for (nc::net::GraphNodeIndex src : graph->AllNodes()) {
      for (nc::net::GraphNodeIndex dst : graph->AllNodes()) {
        if (src == dst) {
          continue;
        }

        double distance_ms =
            duration_cast<milliseconds>(sp.GetDistance(src, dst)).count();
        node_distances_ms.emplace_back(distance_ms);
      }
    }
    node_distances_plot.AddData("", node_distances_ms);

    ctemplate::TemplateDictionary dict("plot");
    dict.SetValue("topology_name", name);

    for (const auto& seed_and_ll : seed_map) {
      double seed = seed_and_ll.first;
      const LLMap& ll_map = seed_and_ll.second;
      std::string subdir = nc::StrCat(root, "/tm_seed_", seed);
      PlotLLMap(name, graph, seed, ll_map, subdir);

      ctemplate::TemplateDictionary* subdict =
          dict.AddSectionDictionary("traffic_matrices");
      subdict->SetValue("tm_seed", nc::StrCat(seed));
    }

    nc::net::GraphStats stats = graph->Stats();
    dict.SetValue("topology_summary", Indent(stats.ToString()));

    PlotAndAddToTemplate(root, "node_distances", node_distances_plot, &dict);

    std::string output;
    ctemplate::ExpandTemplate(kSingleTopTemplate, ctemplate::DO_NOT_STRIP,
                              &dict, &output);
    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root, "/main.rst"));
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

  void PlotDataByLoadAndLocality(
      double load, double locality,
      const std::map<std::string, std::vector<const RoutingConfiguration*>>&
          opt_map,
      const std::string& root) {
    PlotPack plots;
    ctemplate::TemplateDictionary dict("plot");
    std::set<const nc::lp::DemandMatrix*> all_demands;

    for (const auto& optimizer_and_rcs : opt_map) {
      const std::string& opt = optimizer_and_rcs.first;
      const std::vector<const RoutingConfiguration*>& rcs =
          optimizer_and_rcs.second;
      for (const RoutingConfiguration* state : rcs) {
        const TMKey& key = nc::FindOrDieNoPrint(key_by_state_, state);
        const nc::lp::DemandMatrix* demand_matrix =
            nc::FindOrDieNoPrint(demands_by_key_, key);

        all_demands.emplace(demand_matrix);
      }

      HandleSingleOptimizerResult result =
          HandleSingleOptimizer(opt, rcs, &plots);
      ctemplate::TemplateDictionary* subdict =
          dict.AddSectionDictionary("max_ratio_table_rows");
      subdict->SetValue("optimizer", opt);
      subdict->SetValue("value", nc::StrCat(result.max_value));
      subdict->SetValue("link", std::to_string(interesting_states.size()));
      interesting_states.emplace_back(nc::FindOrDie(key_by_state_, result.max));

      subdict = dict.AddSectionDictionary("med_ratio_table_rows");
      subdict->SetValue("optimizer", opt);
      subdict->SetValue("value", nc::StrCat(result.med_value));
      subdict->SetValue("link", std::to_string(interesting_states.size()));
      interesting_states.emplace_back(nc::FindOrDie(key_by_state_, result.med));
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

    PlotAndAddToTemplate(root, "path_ratios", plots.ratios, &dict);
    PlotAndAddToTemplate(root, "path_stretch_rel", plots.path_stretch_rel,
                         &dict);
    PlotAndAddToTemplate(root, "max_path_stretch_rel",
                         plots.path_stretch_max_rel, &dict);
    PlotAndAddToTemplate(root, "path_count", plots.aggregate_path_count, &dict);
    PlotAndAddToTemplate(root, "link_utilization", plots.link_utilization,
                         &dict);
    PlotAndAddToTemplate(root, "link_scales", link_scales_plot, &dict);

    dict.SetValue("tm_load", nc::StrCat(load));
    dict.SetValue("tm_locality", nc::StrCat(locality));
    std::string output;
    ctemplate::ExpandTemplate(kSummaryTemplate, ctemplate::DO_NOT_STRIP, &dict,
                              &output);
    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root, "/main.rst"));
  }

  // Stores graphs.
  std::map<std::string, std::unique_ptr<GraphAndNodeOrder>> graphs_;
  std::map<const nc::net::GraphStorage*, std::string> graph_to_name_;

  // Stores demands.
  std::map<TopologyAndTM, std::unique_ptr<nc::lp::DemandMatrix>>
      demands_by_topology_;
  std::map<TMKey, const nc::lp::DemandMatrix*> demands_by_key_;
  std::map<const RoutingConfiguration*, TMKey> key_by_state_;

  // All data, grouped by load and locality.
  std::map<LoadAndLocality,
           std::map<std::string, std::vector<const RoutingConfiguration*>>>
      data_by_load_and_localiy_;

  // Stores all data.
  TopologyMap top_level_map_;

  // Interesting states to display.
  std::vector<TMKey> interesting_states;
};

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::unique_ptr<DataStorage> data_storage =
      ParseTMGenUtilMetrics(FLAGS_metrics_dir);

  // The data is grouped by optimizer.
  std::vector<std::string> optimizers = nc::Split(FLAGS_optimizers, ",");
  const std::map<std::string, TMStateMap>& data = data_storage->data();

  PlotTotalDelay(data, optimizers);
  PlotCTRRuntime(data);

  // A tree of plots.
  DataPlotter data_plotter;
  for (const std::string& opt : optimizers) {
    const TMStateMap& state_map = nc::FindOrDie(data, opt);
    for (const auto& topology_and_tm_and_rest : state_map) {
      const TopologyAndTM& topology_and_tm = topology_and_tm_and_rest.first;
      const OptimizerTMState* state = topology_and_tm_and_rest.second.get();
      data_plotter.AddData(opt, topology_and_tm, state);
    }
  }

  data_plotter.PlotRoot("plot_tree");
}
