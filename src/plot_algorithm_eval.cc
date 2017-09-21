#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/stats.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/thread_runner.h"
#include "ncode_common/src/viz/grapher.h"
#include "metrics/metrics_parser.h"

DEFINE_string(metrics_dir, "", "The metrics directory.");
DEFINE_string(optimizers, "MinMaxLD,CTRNFC,B4,CTR,MinMaxK10",
              "Optimizers to plot.");
DEFINE_bool(ignore_trivial, true,
            "Ignores topologies/TM combinations where traffic fits on the "
            "shortest path.");
DEFINE_string(focus_stretch_on, "",
              "If empty will plot path stretch for all flows in the data set. "
              "If not will plot only flows in the TM idenified by this.");

static constexpr char kPathStretchMetric[] = "opt_path_stretch_rel";
static constexpr char kPathStretchAbsMetric[] = "opt_path_stretch_ms";
static constexpr char kPathDelayMetric[] = "opt_path_sp_delay_ms";
static constexpr char kPathCountMetric[] = "opt_path_flow_count";
static constexpr char kLinkUtilizationMetric[] = "opt_link_utilization";
static constexpr char kAggregatePathCountMetric[] = "opt_path_count";
static constexpr char kRuntimeMetric[] = "ctr_runtime_ms";
static constexpr char kRuntimeCachedMetric[] = "ctr_runtime_cached_ms";
static constexpr size_t kDiscreteMultiplier = 1000;
static constexpr size_t kPercentilesCount = 10000;

using namespace nc::metrics::parser;

struct OptimizerTMState {
  std::vector<float> path_stretch;
  std::vector<uint32_t> path_stretch_abs_ms;
  std::vector<uint32_t> path_sp_delay_ms;
  std::vector<uint32_t> path_flow_count;
  std::vector<uint32_t> aggregate_path_count;
  std::vector<double> link_utilization;
};

using TopologyAndTM = std::pair<std::string, std::string>;
using TMStateMap = std::map<TopologyAndTM, std::unique_ptr<OptimizerTMState>>;

using DataVector = std::vector<double>;
using DataMap = std::map<std::pair<std::string, std::string>, DataVector>;

class DataStorage {
 public:
  OptimizerTMState* GetTMState(const std::string& topology,
                               const std::string& optimizer,
                               const std::string& tm) {
    std::unique_lock<std::mutex> lock(mu_);
    TMStateMap& tm_state_map = data_[optimizer];
    std::unique_ptr<OptimizerTMState>& tm_state_ptr =
        tm_state_map[{topology, tm}];

    if (!tm_state_ptr) {
      tm_state_ptr = nc::make_unique<OptimizerTMState>();
    }
    return tm_state_ptr.get();
  }

  const std::map<std::string, TMStateMap>& data() const { return data_; }

 private:
  std::map<std::string, TMStateMap> data_;

  // Protects 'data_'.
  std::mutex mu_;
};

class SingleMetricProcessor : public MetricProcessor {
 public:
  using UpdateStateCallback = std::function<void(
      const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state)>;

  SingleMetricProcessor(const std::string& metric, DataStorage* storage,
                        const std::set<TopologyAndTM>* to_ignore,
                        UpdateStateCallback callback)
      : metric_(metric),
        callback_(callback),
        to_ignore_(to_ignore),
        storage_(storage) {}

  bool InterestedInMetric(const std::string& metric_id) override {
    return metric_ == metric_id;
  }

  bool InterestedInFields(const nc::metrics::PBManifestEntry& manifest_entry,
                          uint32_t manifest_index) override {
    // Fields are assumed to be topology,TM,optimizer.
    CHECK(manifest_entry.fields_size() == 3);
    CHECK(manifest_entry.fields(0).type() ==
          nc::metrics::PBMetricField_Type_STRING);
    CHECK(manifest_entry.fields(1).type() ==
          nc::metrics::PBMetricField_Type_STRING);
    CHECK(manifest_entry.fields(2).type() ==
          nc::metrics::PBMetricField_Type_STRING);

    const std::string& topology = manifest_entry.fields(0).string_value();
    const std::string& tm = manifest_entry.fields(1).string_value();
    const std::string& optimizer = manifest_entry.fields(2).string_value();

    if (nc::ContainsKey(*to_ignore_, make_pair(topology, tm))) {
      return false;
    }

    OptimizerTMState* tm_state = storage_->GetTMState(topology, optimizer, tm);

    CHECK(!nc::ContainsKey(manifest_index_to_state_, manifest_index));
    manifest_index_to_state_[manifest_index] = tm_state;
    return true;
  }

  void ProcessEntry(const nc::metrics::PBMetricEntry& entry,
                    const nc::metrics::PBManifestEntry& manifest_entry,
                    uint32_t manifest_index) override {
    nc::Unused(manifest_entry);
    OptimizerTMState* tm_state =
        nc::FindOrDie(manifest_index_to_state_, manifest_index);
    callback_(entry, tm_state);
  }

 private:
  const std::string metric_;

  UpdateStateCallback callback_;

  // Combinations of TM and topology to ignore.
  const std::set<TopologyAndTM>* to_ignore_;

  // The storage that produces OptimizerTMState instances.
  DataStorage* storage_;

  // Maps from a manifest index to state.
  std::map<uint32_t, OptimizerTMState*> manifest_index_to_state_;
};

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

    double total_sp_delay = 0;
    double total_delay = 0;
    for (size_t i = 0; i < tm_state.path_flow_count.size(); ++i) {
      double flow_count = tm_state.path_flow_count[i];
      double abs_stretch = tm_state.path_stretch_abs_ms[i];
      double sp_delay = tm_state.path_sp_delay_ms[i];

      total_sp_delay += sp_delay * flow_count;
      total_delay += (sp_delay + abs_stretch) * flow_count;
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

    double total_delay = 0;
    for (size_t i = 0; i < tm_state.path_flow_count.size(); ++i) {
      double flow_count = tm_state.path_flow_count[i];
      double abs_stretch = tm_state.path_stretch_abs_ms[i];
      double sp_delay = tm_state.path_sp_delay_ms[i];

      total_delay += (sp_delay + abs_stretch) * flow_count;
    }

    out.emplace_back(total_delay);
  }

  return out;
}

// Returns the percentiles of the distribution of a per-path value and the max.
static std::pair<std::vector<double>, std::vector<double>>
GetStretchDistribution(const TMStateMap& tm_state_map) {
  nc::DiscreteDistribution<uint64_t> dist;
  std::vector<double> maxs;

  for (const auto& key_and_data : tm_state_map) {
    const OptimizerTMState& tm_state = *(key_and_data.second);

    double max = 0;
    for (size_t i = 0; i < tm_state.path_flow_count.size(); ++i) {
      double stretch = tm_state.path_stretch[i];
      max = std::max(max, stretch);

      if (FLAGS_focus_stretch_on.empty() ||
          FLAGS_focus_stretch_on == key_and_data.first.second) {
        double abs_stretch = tm_state.path_stretch_abs_ms[i];
        // Have to discretize the value.
        uint64_t value_discrete =
            static_cast<uint64_t>(kDiscreteMultiplier * abs_stretch);
        size_t flow_count = tm_state.path_flow_count[i];
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
    nc::viz::PythonGrapher grapher(out);                           \
    grapher.PlotLine({}, to_plot);                                 \
  }

#define PLOT1D(out, v)                                             \
  {                                                                \
    std::vector<nc::viz::DataSeries1D> to_plot(optimizers.size()); \
    for (size_t i = 0; i < optimizers.size(); ++i) {               \
      to_plot[i] = std::move(v);                                   \
    }                                                              \
    nc::viz::PythonGrapher grapher(out);                           \
    grapher.PlotCDF({}, to_plot);                                  \
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

  nc::viz::PythonGrapher grapher("total_delays");
  grapher.PlotLine({}, to_plot);
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

  std::vector<std::pair<double, double>> xy;
  std::vector<std::pair<double, double>> xy_cached;
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

    xy.emplace_back(size, runtime_data_vector.front());
    xy_cached.emplace_back(size, runtime_data_cached_vector.front());
  }

  std::vector<nc::viz::DataSeries2D> to_plot;
  to_plot.push_back({"Regular", xy_cached});
  to_plot.push_back({"Cold run", xy});
  nc::viz::PythonGrapher grapher("runtime");
  grapher.PlotLine({}, to_plot);
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

  DataStorage data_storage;
  auto path_stretch_processor = nc::make_unique<SingleMetricProcessor>(
      kPathStretchMetric, &data_storage, &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->path_stretch.emplace_back(entry.double_value());
      });

  auto path_stretch_abs_processor = nc::make_unique<SingleMetricProcessor>(
      kPathStretchAbsMetric, &data_storage, &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->path_stretch_abs_ms.emplace_back(entry.uint32_value());
      });

  auto path_delay_processor = nc::make_unique<SingleMetricProcessor>(
      kPathDelayMetric, &data_storage, &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->path_sp_delay_ms.emplace_back(entry.uint32_value());
      });

  auto path_flow_count_processor = nc::make_unique<SingleMetricProcessor>(
      kPathCountMetric, &data_storage, &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->path_flow_count.emplace_back(entry.uint32_value());
      });

  auto link_utilization_processor = nc::make_unique<SingleMetricProcessor>(
      kLinkUtilizationMetric, &data_storage, &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->link_utilization.emplace_back(entry.double_value());
      });

  auto aggregate_path_count_processor = nc::make_unique<SingleMetricProcessor>(
      kAggregatePathCountMetric, &data_storage, &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->aggregate_path_count.emplace_back(entry.uint32_value());
      });

  MetricsParser parser(FLAGS_metrics_dir);
  parser.AddProcessor(std::move(path_stretch_processor));
  parser.AddProcessor(std::move(path_stretch_abs_processor));
  parser.AddProcessor(std::move(path_delay_processor));
  parser.AddProcessor(std::move(path_flow_count_processor));
  parser.AddProcessor(std::move(link_utilization_processor));
  parser.AddProcessor(std::move(aggregate_path_count_processor));
  parser.Parse();

  // The data is grouped by optimizer.
  std::vector<std::string> optimizers = nc::Split(FLAGS_optimizers, ",");
  std::vector<std::unique_ptr<Result>> result_ptrs;
  const std::map<std::string, TMStateMap>& data = data_storage.data();
  for (const std::string& opt : optimizers) {
    LOG(INFO) << "Handle " << opt;
    const TMStateMap& state_map = nc::FindOrDie(data, opt);

    auto result_ptr = HandleSingleOptimizer(opt, state_map);
    result_ptrs.emplace_back(std::move(result_ptr));
  }

  PLOT("path_rel_stretch_out", result_ptrs[i]->path_stretch_rel_data);
  PLOT1D("aggregate_path_count_out", result_ptrs[i]->aggregate_path_count);
  PLOT1D("path_rel_stretch_max_out", result_ptrs[i]->path_stretch_max_rel_data);
  PLOT1D("link_utilization_out", result_ptrs[i]->link_utilization_data);
  PLOT1D("ratios_out", result_ptrs[i]->ratios_data);

  PlotTotalDelay(data, optimizers);
  PlotCTRRuntime(data);
}
