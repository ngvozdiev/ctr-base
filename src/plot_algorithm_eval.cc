#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
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
#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
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
DEFINE_string(
    locality_plot, "",
    "If non-empty will plot locality for the given topology/tm combination.");

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

static void PlotLocality(const TMState& tm_state) {
  std::vector<std::pair<double, double>> values;
  CHECK(tm_state.aggregate_rate_Mbps.size() ==
        tm_state.aggregate_sp_delay_ms.size())
      << tm_state.aggregate_rate_Mbps.size() << " vs "
      << tm_state.aggregate_sp_delay_ms.size();
  for (size_t i = 0; i < tm_state.aggregate_rate_Mbps.size(); ++i) {
    values.emplace_back(tm_state.aggregate_sp_delay_ms[i],
                        tm_state.aggregate_rate_Mbps[i]);
  }

  std::sort(values.begin(), values.end());
  double total = 0;
  for (size_t i = 0; i < values.size(); ++i) {
    std::pair<double, double>& delay_and_rate = values[i];
    total += delay_and_rate.second;
    delay_and_rate.second = total;
  }

  nc::viz::PlotParameters2D params;
  params.title = "SP delay vs cumulative rate";
  params.x_label = "SP delay (ms)";
  params.y_label = "Total cumulative rate (Mbps)";
  nc::viz::LinePlot plot(params);
  plot.AddData("", values);
  plot.PlotToDir(nc::StrCat(FLAGS_output_prefix, "locality_plot"));
}

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

  if (!FLAGS_locality_plot.empty()) {
    std::vector<std::string> split = nc::Split(FLAGS_locality_plot, ":");
    CHECK(split.size() == 2);

    const auto& tm_state_data = data_storage->tm_state_data();
    const TMState& tm_state =
        *(nc::FindOrDieNoPrint(tm_state_data, {split[0], split[1]}));
    PlotLocality(tm_state);
    return 0;
  }

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
}
