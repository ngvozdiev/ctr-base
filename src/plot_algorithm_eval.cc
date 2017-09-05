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

DEFINE_string(metric_dir, "", "The metrics directory.");
DEFINE_string(focus_on, "", "A field to plot absolute values of.");
DEFINE_string(optimizers, "MinMaxLD,CTRNFC,B4,CTR,MinMaxK10",
              "Comma-separated list of optimizers");
DEFINE_bool(ignore_trivial, true,
            "Ignores topologies/TM combinations where traffic fits on the "
            "shortest path.");

static constexpr char kPathStretchRelMetric[] = "opt_path_stretch_rel";
static constexpr char kPathStretchAbsMetric[] = "opt_path_stretch_ms";
static constexpr char kPathDelayMetric[] = "opt_path_sp_delay_ms";
static constexpr char kPathCountMetric[] = "opt_path_flow_count";
static constexpr char kLinkUtilizationMetric[] = "opt_link_utilization";
static constexpr char kAggregatePathCountMetric[] = "opt_path_count";
static constexpr size_t kDiscreteMultiplier = 1000;
static constexpr size_t kPercentilesCount = 10000;

using DataVector = std::vector<double>;
using DataMap = const std::map<std::pair<std::string, std::string>, DataVector>;

using namespace nc::metrics::parser;

struct OptimizerTMState {
  // Identifies the TM.
  std::string id;

  std::vector<float> path_stretch;
  std::vector<uint32_t> path_stretch_abs_ms;
  std::vector<uint32_t> path_sp_delay_ms;
  std::vector<uint32_t> path_flow_count;
  std::vector<uint32_t> aggregate_path_count;
  std::vector<double> link_utilization;
};

class OptimizerParser : public MetricProcessor {
  bool InterestedInFields(const nc::metrics::PBManifestEntry& manifest_entry,
                          uint32_t manifest_index) override {
    manifest_entry.
  }

 private:
  std::string opt;
};

// Returns the set of fields for which all aggregates use only one path.
static std::set<std::string> FieldsWithSinglePath(
    const DataMap& path_count_data) {
  std::set<std::string> out;
  for (const auto& key_and_data : path_count_data) {
    const std::pair<std::string, std::string>& metric_and_fields =
        key_and_data.first;
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
      out.emplace(metric_and_fields.second);
    }
  }

  return out;
}

static std::vector<double> GetPathRatios(const std::set<std::string>& to_ignore,
                                         const DataMap& flow_counts_data,
                                         const DataMap& path_stretch_abs_data,
                                         const DataMap& path_sp_delay_data) {
  std::vector<double> out;
  for (const auto& key_and_data : flow_counts_data) {
    const std::pair<std::string, std::string>& metric_and_fields =
        key_and_data.first;
    const std::string& fields = metric_and_fields.second;
    if (nc::ContainsKey(to_ignore, fields)) {
      continue;
    }

    const DataVector& flow_counts = key_and_data.second;
    const DataVector& path_stretch_abs = nc::FindOrDieNoPrint(
        path_stretch_abs_data, {kPathStretchAbsMetric, fields});
    const DataVector& path_sp_delay =
        nc::FindOrDieNoPrint(path_sp_delay_data, {kPathDelayMetric, fields});
    CHECK(flow_counts.size() == path_stretch_abs.size());
    CHECK(flow_counts.size() == path_sp_delay.size());

    double total_sp_delay = 0;
    double total_delay = 0;
    for (size_t i = 0; i < flow_counts.size(); ++i) {
      double abs_stretch = path_stretch_abs[i];
      double sp_delay = path_sp_delay[i];

      total_sp_delay += sp_delay * flow_counts[i];
      total_delay += (sp_delay + abs_stretch) * flow_counts[i];
    }

    double change = (total_delay - total_sp_delay) / total_sp_delay;
    out.emplace_back(change);
  }

  return out;
}

// Returns the percentiles of the distribution of a per-path value and the max.
static std::pair<std::vector<double>, std::vector<double>> Parse(
    const std::set<std::string>& to_ignore, const DataMap& path_counts_data,
    const DataMap& per_path_data) {
  nc::DiscreteDistribution<uint64_t> dist;
  std::vector<double> maxs;

  for (const auto& key_and_data : per_path_data) {
    const std::pair<std::string, std::string>& metric_and_fields =
        key_and_data.first;
    if (nc::ContainsKey(to_ignore, metric_and_fields.second)) {
      continue;
    }

    const DataVector& path_data = key_and_data.second;

    // Need to get the flow counts for each path from the kPathCountMetric. The
    // fields will be the same.
    const DataVector* flow_counts = nc::FindOrNull(
        path_counts_data, {kPathCountMetric, metric_and_fields.second});
    if (flow_counts == nullptr) {
      LOG(ERROR) << "No flow counts for " << metric_and_fields.second;
      continue;
    }

    // The two data vectors should have the same number of elements---one for
    // each path with either the data or the flow count on the path.
    CHECK(flow_counts->size() == path_data.size())
        << metric_and_fields.second << " " << flow_counts->size() << " vs "
        << path_data.size();

    double max = 0;
    for (size_t i = 0; i < flow_counts->size(); ++i) {
      // Have to discretize the value.
      double value = path_data[i];
      max = std::max(max, value);

      uint64_t value_discrete =
          static_cast<uint64_t>(kDiscreteMultiplier * value);

      // The flow count should always be integer.
      double flow_count;
      CHECK(std::modf((*flow_counts)[i], &flow_count) == 0.0);
      dist.Add(value_discrete, static_cast<size_t>(flow_count));
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

static std::vector<double> ParseSingle(const DataVector& path_counts,
                                       const DataVector& path_data,
                                       double* total, double* count) {
  CHECK(path_counts.size() == path_data.size());
  nc::DiscreteDistribution<uint64_t> dist;
  for (size_t i = 0; i < path_counts.size(); ++i) {
    double value = path_data[i];
    uint64_t value_discrete =
        static_cast<uint64_t>(kDiscreteMultiplier * value);

    double path_count;
    CHECK(std::modf(path_counts[i], &path_count) == 0.0);
    dist.Add(value_discrete, static_cast<size_t>(path_count));
  }

  std::vector<uint64_t> percentiles = dist.Percentiles(kPercentilesCount);
  CHECK(percentiles.size() == kPercentilesCount + 1);

  std::vector<double> out;
  out.reserve(percentiles.size());
  for (uint64_t discrete_value : percentiles) {
    out.emplace_back(static_cast<double>(discrete_value) / kDiscreteMultiplier);
  }

  if (total != nullptr) {
    *total = dist.summary_stats().sum() / kDiscreteMultiplier;
  }

  if (count != nullptr) {
    *count = dist.summary_stats().count();
  }

  return out;
}

struct Result {
  nc::viz::DataSeries2D path_stretch_rel_data;
  nc::viz::DataSeries1D path_stretch_max_rel_data;
  nc::viz::DataSeries1D link_utilization_data;
  nc::viz::DataSeries1D ratios_data;
  nc::viz::DataSeries2D aggregate_path_count;
};

static nc::viz::DataSeries1D FlattenData(const std::set<std::string>& to_ignore,
                                         const DataMap& data_map,
                                         const std::string& label) {
  std::vector<double> values;
  for (const auto& key_and_value : data_map) {
    const std::pair<std::string, std::string>& metric_and_fields =
        key_and_value.first;
    if (nc::ContainsKey(to_ignore, metric_and_fields.second)) {
      continue;
    }

    for (const auto& v : key_and_value.second) {
      values.emplace_back(v);
    }
  }

  nc::viz::DataSeries1D data_series;
  data_series.data = std::move(values);
  data_series.label = label;

  return data_series;
}

static std::unique_ptr<Result> HandleSingleOptimizer(
    const std::set<std::string>& to_ignore, const std::string& optimizer) {
  std::string fields_regex = nc::StrCat(".*", optimizer, "$");

  LOG(INFO) << "Parsing metric " << kPathCountMetric << " fields "
            << fields_regex;
  DataMap flow_counts = SimpleParseNumericDataNoTimestamps(
      FLAGS_metric_dir, kPathCountMetric, fields_regex);
  LOG(INFO) << "Done";

  auto result_ptr = nc::make_unique<Result>();
  DataMap path_stretch_rel_data = SimpleParseNumericDataNoTimestamps(
      FLAGS_metric_dir, kPathStretchRelMetric, fields_regex);
  DataMap path_stretch_abs_data = SimpleParseNumericDataNoTimestamps(
      FLAGS_metric_dir, kPathStretchAbsMetric, fields_regex);
  DataMap path_sp_delay_data = SimpleParseNumericDataNoTimestamps(
      FLAGS_metric_dir, kPathDelayMetric, fields_regex);

  std::vector<double> ratios = GetPathRatios(
      to_ignore, flow_counts, path_stretch_abs_data, path_sp_delay_data);
  result_ptr->ratios_data.label = optimizer;
  result_ptr->ratios_data.data = std::move(ratios);

  // Will plot the percentiles to get a CDF.
  std::vector<double> percentiles;
  std::vector<double> maxs;
  std::tie(percentiles, maxs) =
      Parse(to_ignore, flow_counts, path_stretch_rel_data);
  for (size_t i = 0; i < percentiles.size(); ++i) {
    double p = percentiles[i];
    result_ptr->path_stretch_rel_data.data.emplace_back(p, i);
  }
  result_ptr->path_stretch_max_rel_data.data = std::move(maxs);
  result_ptr->path_stretch_rel_data.label = optimizer;
  result_ptr->path_stretch_max_rel_data.label = optimizer;

  std::vector<double> per_aggregate_path_counts;
  DataMap per_aggregate_data = SimpleParseNumericDataNoTimestamps(
      FLAGS_metric_dir, kAggregatePathCountMetric, fields_regex);
  for (const auto& key_and_data : per_aggregate_data) {
    const std::pair<std::string, std::string>& metric_and_fields =
        key_and_data.first;
    if (nc::ContainsKey(to_ignore, metric_and_fields.second)) {
      continue;
    }

    const DataVector& path_data = key_and_data.second;
    for (const auto& datapoint : path_data) {
      per_aggregate_path_counts.emplace_back(datapoint);
    }
  }
  percentiles = nc::Percentiles(&per_aggregate_path_counts, kPercentilesCount);
  for (size_t i = 0; i < percentiles.size(); ++i) {
    double p = percentiles[i];
    result_ptr->aggregate_path_count.data.emplace_back(p, i);
  }
  result_ptr->aggregate_path_count.label = optimizer;

  DataMap link_utilization_data = SimpleParseNumericDataNoTimestamps(
      FLAGS_metric_dir, kLinkUtilizationMetric, fields_regex);
  result_ptr->link_utilization_data =
      FlattenData(to_ignore, link_utilization_data, optimizer);

  return result_ptr;
}

static std::unique_ptr<nc::viz::DataSeries2D> HandleFocusOn(
    const std::string& optimizer) {
  std::string field_regex =
      nc::StrCat(".*", FLAGS_focus_on, ".*", optimizer, "$");
  LOG(INFO) << "Parsing " << kPathCountMetric << " for " << field_regex;
  DataMap path_count_data = SimpleParseNumericDataNoTimestamps(
      FLAGS_metric_dir, kPathCountMetric, field_regex);

  LOG(INFO) << "Parsing " << kPathStretchAbsMetric << " for " << field_regex;
  DataMap path_stretch_abs_data = SimpleParseNumericDataNoTimestamps(
      FLAGS_metric_dir, kPathStretchAbsMetric, field_regex);

  // The two datamaps should only have one value---the fields that we are
  // focusing on.
  CHECK(path_count_data.size() == 1);
  CHECK(path_stretch_abs_data.size() == 1);

  const DataVector& path_count_data_vector = path_count_data.begin()->second;
  const DataVector& path_stretch_data_vector =
      path_stretch_abs_data.begin()->second;

  double total;
  double count;
  std::vector<double> percentiles = ParseSingle(
      path_count_data_vector, path_stretch_data_vector, &total, &count);
  LOG(INFO) << "For " << optimizer << " total " << total << " count " << count;

  auto to_plot = nc::make_unique<nc::viz::DataSeries2D>();
  to_plot->label = optimizer;
  for (size_t i = 0; i < percentiles.size(); ++i) {
    double p = percentiles[i];
    to_plot->data.emplace_back(p, i);
  }

  return to_plot;
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

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_metric_dir.empty()) << "need --metric_file";

  std::vector<std::string> optimizers = nc::Split(FLAGS_optimizers, ",");
  if (!FLAGS_focus_on.empty()) {
    std::vector<std::unique_ptr<nc::viz::DataSeries2D>> result_ptrs =
        nc::RunInParallelWithResult<std::string, nc::viz::DataSeries2D>(
            optimizers, [](const std::string& opt) {
              return HandleFocusOn(opt);
            }, optimizers.size());

    PLOT("path_abs_stretch_out_focus", *result_ptrs[i]);
    return 0;
  }

  std::set<std::string> to_ignore;
  if (FLAGS_ignore_trivial) {
    std::string fields_regex = nc::StrCat(".*B4$");
    DataMap path_count_data = SimpleParseNumericDataNoTimestamps(
        FLAGS_metric_dir, kAggregatePathCountMetric, fields_regex);
    to_ignore = FieldsWithSinglePath(path_count_data);
    LOG(INFO) << "Will ignore " << to_ignore.size()
              << " cases satisfiable with one path";
  }

  std::vector<std::unique_ptr<Result>> result_ptrs =
      nc::RunInParallelWithResult<std::string, Result>(
          optimizers, [&to_ignore](const std::string& opt) {
            return HandleSingleOptimizer(to_ignore, opt);
          }, 1);

  PLOT("path_rel_stretch_out", result_ptrs[i]->path_stretch_rel_data);
  PLOT("aggregate_path_count_out", result_ptrs[i]->aggregate_path_count);
  PLOT1D("path_rel_stretch_max_out", result_ptrs[i]->path_stretch_max_rel_data);
  PLOT1D("link_utilization_out", result_ptrs[i]->link_utilization_data);
  PLOT1D("ratios_out", result_ptrs[i]->ratios_data);

  DataMap runtime_data = SimpleParseNumericDataNoTimestamps(
      FLAGS_metric_dir, "ctr_runtime_ms", ".*");
  DataMap runtime_data_cached = SimpleParseNumericDataNoTimestamps(
      FLAGS_metric_dir, "ctr_runtime_cached_ms", ".*");
  nc::viz::PythonGrapher grapher("runtime_out");
  grapher.PlotCDF(
      {}, {FlattenData(to_ignore, runtime_data, "Runtime (ms)"),
           FlattenData(to_ignore, runtime_data_cached, "Runtime cached (ms)")});
}
