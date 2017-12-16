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

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/algorithm.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/viz/grapher.h"
#include "metrics/metrics_parser.h"
#include "plot_algorithm_eval_tools.h"

DEFINE_string(metrics_dir, "", "The metrics directory.");
DEFINE_string(optimizers, "MinMaxLD,CTRNFC,B4,CTR,MinMaxK10",
              "Optimizers to plot.");
DEFINE_bool(ignore_trivial, true,
            "Ignores topologies/TM combinations where traffic fits on the "
            "shortest path.");
DEFINE_string(output_prefix, "",
              "Prefix that will be added to all output directories.");

static constexpr char kAggregatePathCountMetric[] = "opt_path_count";
static constexpr char kRuntimeMetric[] = "ctr_runtime_ms";
static constexpr char kRuntimeCachedMetric[] = "ctr_runtime_cached_ms";

using namespace nc::metrics::parser;
using namespace ctr::alg_eval;

// Returns the set of fields for which all aggregates use only one path.
static std::set<TopologyAndTM> FieldsWithSinglePath(
    const DataMap& path_count_data) {
  std::set<TopologyAndTM> out;
  for (const auto& key_and_data : path_count_data) {
    const std::pair<std::string, std::string>& metric_and_fields =
        key_and_data.first;

    // Assuming the field string is <topology_file>:<tm_file>:<optimizer>
    std::string field_string = metric_and_fields.second;
    std::vector<std::string> split = nc::Split(field_string, ":");
    CHECK(split.size() == 3);

    const DataVector& path_counts = key_and_data.second;
    bool all_ones = true;
    for (double count : path_counts) {
      double path_count;
      CHECK(std::modf(count, &path_count) == 0.0);

      if (path_count != 1.0) {
        all_ones = false;
        break;
      }
    }

    if (all_ones) {
      out.emplace(split[0], split[1]);
    }
  }

  return out;
}

static std::vector<double> GetPathRatios(const TMStateMap& tm_state_map,
                                         const std::string* topology) {
  std::vector<double> out;
  for (const auto& key_and_data : tm_state_map) {
    const OptimizerTMState& tm_state = *(key_and_data.second);

    if (topology != nullptr) {
      const TopologyAndTM& topology_and_tm = key_and_data.first;
      if (*topology != topology_and_tm.first) {
        continue;
      }
    }

    std::vector<AggregateTMState> aggregates = tm_state.GetAggregates();
    double total_sp_delay = 0;
    double total_delay = 0;
    for (const auto& aggregate : aggregates) {
      double sp_delay_ms = aggregate.sp_delay_ms;
      for (const auto& delay_and_count : aggregate.paths) {
        double path_delay_ms = delay_and_count.first;
        double abs_stretch = path_delay_ms - sp_delay_ms;
        uint32_t flow_count = delay_and_count.second;
        total_sp_delay += sp_delay_ms * flow_count;
        total_delay += (sp_delay_ms + abs_stretch) * flow_count;
      }
    }

    double change = (total_delay - total_sp_delay) / total_sp_delay;
    out.emplace_back(change);
  }

  return out;
}

static std::pair<double, std::string> MedianPathRatio(
    const std::vector<double>& ratios, const TMStateMap& tm_state_map) {
  std::vector<std::pair<double, std::string>> to_sort;
  size_t i = 0;
  for (const auto& key_and_data : tm_state_map) {
    to_sort.emplace_back(ratios[i], key_and_data.first.second);

    ++i;
  }

  std::sort(to_sort.begin(), to_sort.end());
  return to_sort[to_sort.size() / 2];
}

static std::vector<double> GetTotalDelays(const TMStateMap& tm_state_map) {
  std::vector<double> out;
  for (const auto& key_and_data : tm_state_map) {
    const OptimizerTMState& tm_state = *(key_and_data.second);
    std::vector<AggregateTMState> aggregates = tm_state.GetAggregates();
    double total_delay = 0;
    for (const auto& aggregate : aggregates) {
      double sp_delay_ms = aggregate.sp_delay_ms;
      for (const auto& delay_and_count : aggregate.paths) {
        double path_delay_ms = delay_and_count.first;
        double abs_stretch = path_delay_ms - sp_delay_ms;
        uint32_t flow_count = delay_and_count.second;
        total_delay += (sp_delay_ms + abs_stretch) * flow_count;
      }
    }

    out.emplace_back(total_delay);
  }

  return out;
}

struct Result {
  nc::viz::DataSeries2D path_stretch_rel_data;
  nc::viz::DataSeries1D path_stretch_max_rel_data;
  nc::viz::DataSeries1D link_utilization_data;
  nc::viz::DataSeries1D ratios_data;
  nc::viz::DataSeries1D aggregate_path_count;
};

static std::unique_ptr<Result> HandleSingleOptimizer(
    const std::string& optimizer, const TMStateMap& tm_state_map) {
  auto result_ptr = nc::make_unique<Result>();
  std::vector<double> ratios = GetPathRatios(tm_state_map, nullptr);
  std::pair<double, std::string> median_ratio =
      MedianPathRatio(ratios, tm_state_map);
  LOG(INFO) << "Opt " << optimizer << " median " << median_ratio.first
            << " top " << median_ratio.second;

  result_ptr->ratios_data.label = optimizer;
  result_ptr->ratios_data.data = std::move(ratios);

  // Will plot the percentiles to get a CDF.
  std::vector<double> percentiles;
  std::vector<double> maxs;
  std::tie(percentiles, maxs) = GetStretchDistribution(tm_state_map);
  for (size_t i = 0; i < percentiles.size(); ++i) {
    double p = percentiles[i];
    result_ptr->path_stretch_rel_data.data.emplace_back(p, i);
  }
  result_ptr->path_stretch_max_rel_data.data = std::move(maxs);
  result_ptr->path_stretch_rel_data.label = optimizer;
  result_ptr->path_stretch_max_rel_data.label = optimizer;

  std::vector<double> per_aggregate_path_counts;
  for (const auto& key_and_data : tm_state_map) {
    const OptimizerTMState& tm_state = *(key_and_data.second);
    for (uint32_t path_count : tm_state.aggregate_path_count) {
      per_aggregate_path_counts.emplace_back(path_count);
    }
  }
  result_ptr->aggregate_path_count.data = std::move(per_aggregate_path_counts);
  result_ptr->aggregate_path_count.label = optimizer;

  std::vector<double> link_utilizations;
  for (const auto& key_and_data : tm_state_map) {
    const OptimizerTMState& tm_state = *(key_and_data.second);
    for (float link_utilization : tm_state.link_utilization) {
      link_utilizations.emplace_back(link_utilization);
    }
  }
  result_ptr->link_utilization_data.data = std::move(link_utilizations);
  result_ptr->link_utilization_data.label = optimizer;

  return result_ptr;
}

#define PLOT(out, v)                                               \
  {                                                                \
    std::vector<nc::viz::DataSeries2D> to_plot(optimizers.size()); \
    for (size_t i = 0; i < optimizers.size(); ++i) {               \
      to_plot[i] = std::move(v);                                   \
    }                                                              \
    nc::viz::LinePlot plot;                                        \
    plot.AddData(to_plot);                                         \
    plot.PlotToDir(out);                                           \
  }

#define PLOT1D(out, v)                                             \
  {                                                                \
    std::vector<nc::viz::DataSeries1D> to_plot(optimizers.size()); \
    for (size_t i = 0; i < optimizers.size(); ++i) {               \
      to_plot[i] = std::move(v);                                   \
    }                                                              \
    nc::viz::CDFPlot plot;                                         \
    plot.AddData(to_plot);                                         \
    plot.PlotToDir(out);                                           \
  }

static void PlotTotalDelay(const std::map<std::string, TMStateMap>& data,
                           const std::vector<std::string>& optimizers) {
  std::vector<nc::viz::DataSeries2D> to_plot;
  for (const std::string& opt : optimizers) {
    const TMStateMap& state_map = nc::FindOrDie(data, opt);

    std::vector<std::pair<double, double>> xy;
    std::vector<double> total_delays = GetTotalDelays(state_map);
    for (size_t i = 0; i < total_delays.size(); ++i) {
      xy.emplace_back(i, total_delays[i]);
    }

    to_plot.push_back({opt, xy});
  }

  nc::viz::LinePlot plot;
  plot.AddData(to_plot);
  plot.PlotToDir("total_delays");
}

// Returns a list of topology/tm pairs ordered by number of aggregates.
static std::vector<std::pair<size_t, TopologyAndTM>>
TopologiesAndTMOrderedBySize(const TMStateMap& data) {
  std::vector<std::pair<size_t, TopologyAndTM>> out;
  for (const auto& topology_and_tm_and_state : data) {
    const TopologyAndTM& topology_and_tm = topology_and_tm_and_state.first;
    const OptimizerTMState& tm_state = *(topology_and_tm_and_state.second);

    size_t aggregate_count = tm_state.aggregate_path_count.size();
    out.emplace_back(aggregate_count, topology_and_tm);
  }

  std::sort(out.begin(), out.end());
  return out;
}

static void PlotCTRRuntime(const std::map<std::string, TMStateMap>& data) {
  DataMap runtime_data = SimpleParseNumericDataNoTimestamps(
      FLAGS_metrics_dir, kRuntimeMetric, ".*");
  DataMap runtime_data_cached = SimpleParseNumericDataNoTimestamps(
      FLAGS_metrics_dir, kRuntimeCachedMetric, ".*");

  const TMStateMap& state_map = nc::FindOrDie(data, "CTR");
  std::vector<std::pair<size_t, TopologyAndTM>> runs_ordered =
      TopologiesAndTMOrderedBySize(state_map);

  std::vector<double> y;
  std::vector<double> y_cached;
  for (const auto& size_and_topology : runs_ordered) {
    size_t size = size_and_topology.first;
    const TopologyAndTM& topology_and_tm = size_and_topology.second;

    std::string key =
        nc::StrCat(topology_and_tm.first, ":", topology_and_tm.second);
    const DataVector& runtime_data_vector =
        nc::FindOrDieNoPrint(runtime_data, std::make_pair(kRuntimeMetric, key));
    const DataVector& runtime_data_cached_vector = nc::FindOrDieNoPrint(
        runtime_data_cached, std::make_pair(kRuntimeCachedMetric, key));

    CHECK(runtime_data_vector.size() == 1);
    CHECK(runtime_data_cached_vector.size() == 1);

    y.emplace_back(runtime_data_vector.front());
    y_cached.emplace_back(runtime_data_cached_vector.front());
  }

  nc::viz::CDFPlot plot;
  plot.AddData("Regular", y_cached);
  plot.AddData("Cold run", y);
  plot.PlotToDir("runtime");
}

class DataPlotter {
 public:
  void AddData(const std::string& opt, const TopologyAndTM& topology_and_tm,
               const OptimizerTMState* optimizer_state) {
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

    double seed;
    std::string seed_str = nc::FindOrDie(demand_matrix->properties(), "seed");
    CHECK(nc::safe_strtod(seed_str, &seed));

    TMKey key = std::make_tuple(graph, seed, locality);
    demands_by_key_[key] = demand_matrix.get();
    top_level_map_[graph][seed][locality][opt] = optimizer_state;
  }

  void PlotRoot(const std::string& root) {
    for (const auto& graph_and_seed_map : top_level_map_) {
      const nc::net::GraphStorage* graph = graph_and_seed_map.first;
      const SeedMap& seed_map = graph_and_seed_map.second;
      const std::string& name = nc::FindOrDie(graph_to_name_, graph);

      std::string subdir = nc::StrCat(root, "/", name);
      PlotSeedMap(graph, seed_map, subdir);
    }
  }

 private:
  // Graph, seed and locality.
  using TMKey = std::tuple<const nc::net::GraphStorage*, double, double>;

  struct GraphAndNodeOrder {
    GraphAndNodeOrder(const nc::net::GraphBuilder& builder,
                      const std::vector<std::string>& node_order)
        : graph(builder), node_order(node_order) {}

    nc::net::GraphStorage graph;
    std::vector<std::string> node_order;
  };

  // Maps an optimizer to its state.
  using OptMap = std::map<std::string, const OptimizerTMState*>;

  // Maps locality to all traffic matrices with the same locality.
  using LocalityMap = std::map<double, OptMap>;

  // Groups TMs with the same seed value.
  using SeedMap = std::map<double, LocalityMap>;

  // Maps a topology to all traffic matrices on the same topology.
  using TopologyMap = std::map<const nc::net::GraphStorage*, SeedMap>;

  // Plots different optimizers for the same TM.
  void PlotOptMap(const OptMap& opt_map, const std::string& root) {
    nc::viz::LinePlot abs_path_stretch_plot(
        {"Absolute path stretch", "milliseconds longer than SP", "CDF"});
    nc::viz::LinePlot rel_path_stretch_plot(
        {"Relative path stretch", "fraction change over SP", "CDF"});
    nc::viz::CDFPlot path_count_plot(
        {"Number of paths per aggregate", "path count"});
    nc::viz::CDFPlot link_utilization_plot({"Link utilization", "utilization"});

    for (const auto& optimizer_and_state : opt_map) {
      const std::string& opt = optimizer_and_state.first;
      const OptimizerTMState* state = optimizer_and_state.second;
      CHECK(state != nullptr);

      std::vector<double> abs_percentiles =
          GetStretchDistribution(*state, true);
      std::vector<std::pair<double, double>> abs_to_plot;
      for (size_t i = 0; i < abs_percentiles.size(); ++i) {
        abs_to_plot.emplace_back(abs_percentiles[i], i);
      }

      std::vector<double> rel_percentiles =
          GetStretchDistribution(*state, false);
      std::vector<std::pair<double, double>> rel_to_plot;
      for (size_t i = 0; i < rel_percentiles.size(); ++i) {
        rel_to_plot.emplace_back(rel_percentiles[i], i);
      }

      abs_path_stretch_plot.AddData(opt, abs_to_plot);
      rel_path_stretch_plot.AddData(opt, rel_to_plot);
      path_count_plot.AddData(opt, state->aggregate_path_count);
      link_utilization_plot.AddData(opt, state->link_utilization);
    }

    abs_path_stretch_plot.PlotToDir(nc::StrCat(root, "/abs_path_stretch"));
    rel_path_stretch_plot.PlotToDir(nc::StrCat(root, "/rel_path_stretch"));
    path_count_plot.PlotToDir(nc::StrCat(root, "/path_count"));
    link_utilization_plot.PlotToDir(nc::StrCat(root, "/link_utilization"));
  }

  void PlotLocalityMap(const nc::net::GraphStorage* graph, double seed,
                       const LocalityMap& tm_map, const std::string& root) {
    using namespace std::chrono;

    nc::viz::CDFPlot demand_sizes_plot(
        {"Demand sizes (all possible demands considered)", "size (Mbps)"});
    nc::viz::CDFPlot sp_utilizations_plot(
        {"Shortest path link utilization", "link utilization"});
    nc::viz::LinePlot cumulative_demands_plot(
        {"Cumulative distance vs aggregate size", "SP distance (ms)",
         "Cumulative size (Mbps)"});

    for (const auto& locality_and_opt_map : tm_map) {
      double locality = locality_and_opt_map.first;
      const OptMap& opt_map = locality_and_opt_map.second;
      const nc::lp::DemandMatrix* demand_matrix = nc::FindOrDieNoPrint(
          demands_by_key_, std::make_tuple(graph, seed, locality));

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
      CHECK(demand_sizes_Mbps.size() < possible_count);
      std::string id = nc::StrCat("locality ", locality);

      demand_sizes_plot.AddData(id, demand_sizes_Mbps);
      sp_utilizations_plot.AddData(id, sp_utilizations);
      cumulative_demands_plot.AddData(id, cumulative_distances);

      // Need to plot the per-optimizer data as well. Will do so in a separate
      // subdir for each TM.
      std::string subdir = nc::StrCat(root, "/tm_locality_", locality);
      PlotOptMap(opt_map, subdir);
    }

    demand_sizes_plot.PlotToDir(nc::StrCat(root, "/demand_sizes"));
    sp_utilizations_plot.PlotToDir(nc::StrCat(root, "/sp_utilization"));
    cumulative_demands_plot.PlotToDir(nc::StrCat(root, "/cumulative_demands"));
  }

  void PlotSeedMap(const nc::net::GraphStorage* graph, const SeedMap& seed_map,
                   const std::string& root) {
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
    node_distances_plot.PlotToDir(nc::StrCat(root, "/node_distances"));

    for (const auto& seed_and_locality : seed_map) {
      double seed = seed_and_locality.first;
      const LocalityMap& locality_map = seed_and_locality.second;

      std::string subdir = nc::StrCat(root, "/tm_seed_", seed);
      PlotLocalityMap(graph, seed, locality_map, subdir);
    }
  }

  // Stores graphs.
  std::map<std::string, std::unique_ptr<GraphAndNodeOrder>> graphs_;
  std::map<const nc::net::GraphStorage*, std::string> graph_to_name_;

  // Stores demands.
  std::map<TopologyAndTM, std::unique_ptr<nc::lp::DemandMatrix>>
      demands_by_topology_;
  std::map<TMKey, const nc::lp::DemandMatrix*> demands_by_key_;

  // Stores all data.
  TopologyMap top_level_map_;
};

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_metrics_dir.empty()) << "need --metrics_dir";

  std::set<TopologyAndTM> to_ignore;
  if (FLAGS_ignore_trivial) {
    std::string fields_regex = nc::StrCat(".*B4$");
    DataMap path_count_data = SimpleParseNumericDataNoTimestamps(
        FLAGS_metrics_dir, kAggregatePathCountMetric, fields_regex);
    to_ignore = FieldsWithSinglePath(path_count_data);
    LOG(INFO) << "Will ignore " << to_ignore.size()
              << " cases satisfiable with one path";
  }

  std::unique_ptr<DataStorage> data_storage =
      ParseTMGenUtilMetrics(FLAGS_metrics_dir);

  // The data is grouped by optimizer.
  std::vector<std::string> optimizers = nc::Split(FLAGS_optimizers, ",");
  std::vector<std::unique_ptr<Result>> result_ptrs;
  const std::map<std::string, TMStateMap>& data = data_storage->data();
  for (const std::string& opt : optimizers) {
    LOG(INFO) << "Handle " << opt;
    const TMStateMap& state_map = nc::FindOrDie(data, opt);

    std::set<std::string> topologies;
    for (const auto& topolog_and_tm_and_rest : state_map) {
      const std::string topology = topolog_and_tm_and_rest.first.first;
      topologies.emplace(topology);
    }
    LOG(INFO) << "Top count " << topologies.size();

    auto result_ptr = HandleSingleOptimizer(opt, state_map);
    result_ptrs.emplace_back(std::move(result_ptr));
  }

  PLOT(nc::StrCat(FLAGS_output_prefix, "path_rel_stretch_out"),
       result_ptrs[i]->path_stretch_rel_data);
  PLOT1D(nc::StrCat(FLAGS_output_prefix, "aggregate_path_count_out"),
         result_ptrs[i]->aggregate_path_count);
  PLOT1D(nc::StrCat(FLAGS_output_prefix, "path_rel_stretch_max_out"),
         result_ptrs[i]->path_stretch_max_rel_data);
  PLOT1D(nc::StrCat(FLAGS_output_prefix, "link_utilization_out"),
         result_ptrs[i]->link_utilization_data);
  PLOT1D(nc::StrCat(FLAGS_output_prefix, "ratios_out"),
         result_ptrs[i]->ratios_data);

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
