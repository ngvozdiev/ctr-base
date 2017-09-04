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

DEFINE_string(metric_file, "", "The metrics file.");
DEFINE_string(focus_on, "", "A field to plot absolute values of.");
DEFINE_string(optimizers, "MinMax,MinMaxLD,CTRNFC,B4,CTR,MinMaxK10",
              "Comma-separated list of optimizers");

static constexpr char kPathStretchRelMetric[] = "opt_path_stretch_rel";
static constexpr char kPathStretchAbsMetric[] = "opt_path_stretch_ms";
static constexpr char kPathCountMetric[] = "opt_path_flow_count";
static constexpr char kPathPenaltyMetric[] = "opt_path_unmet_demand_bps";
static constexpr char kAggregatePathCountMetric[] = "opt_path_count";
static constexpr size_t kDiscreteMultiplier = 1000;
static constexpr size_t kPercentilesCount = 10000;
static constexpr uint64_t kMaxUint = std::numeric_limits<uint64_t>::max();

using DataVector = std::vector<std::pair<uint64_t, double>>;
using DataMap = const std::map<std::pair<std::string, std::string>, DataVector>;

using namespace nc::metrics::parser;

// Returns the percentiles of the distribution of a per-path value.
static std::vector<double> Parse(const DataMap& path_counts_data,
                                 const DataMap& per_path_data,
                                 const DataMap& per_path_penalty,
                                 double penalty) {
  nc::DiscreteDistribution<uint64_t> dist;

  for (const auto& key_and_data : per_path_data) {
    const std::pair<std::string, std::string>& metric_and_fields =
        key_and_data.first;
    const DataVector& path_data = key_and_data.second;

    // Need to get the flow counts for each path from the kPathCountMetric. The
    // fields will be the same.
    const DataVector* path_counts = nc::FindOrNull(
        path_counts_data, {kPathCountMetric, metric_and_fields.second});
    if (path_counts == nullptr) {
      LOG(ERROR) << "No path counts for " << metric_and_fields.second;
      continue;
    }

    // Path penalties, a fixed penalty will be applied for all values for which
    // the penalty is > 1.
    const DataVector* path_penalties = nc::FindOrNull(
        per_path_penalty, {kPathPenaltyMetric, metric_and_fields.second});

    // The two data vectors should have the same number of elements---one for
    // each path with either the data or the flow count on the path.
    CHECK(path_counts->size() == path_data.size())
        << metric_and_fields.second << " " << path_counts->size() << " vs "
        << path_data.size();
    if (path_penalties != nullptr) {
      CHECK(path_penalties->size() == path_data.size());
    }

    for (size_t i = 0; i < path_counts->size(); ++i) {
      double extra_penalty = 0;
      if (path_penalties != nullptr) {
        double p = (*path_penalties)[i].second;
        if (p > 1.0) {
          extra_penalty += penalty;
        }
      }

      // Have to discretize the value.
      double value = path_data[i].second + extra_penalty;
      uint64_t value_discrete =
          static_cast<uint64_t>(kDiscreteMultiplier * value);

      // The path count should always be integer.
      double path_count;
      CHECK(std::modf((*path_counts)[i].second, &path_count) == 0.0);
      dist.Add(value_discrete, static_cast<size_t>(path_count));
    }
  }

  std::vector<uint64_t> percentiles = dist.Percentiles(kPercentilesCount);
  CHECK(percentiles.size() == kPercentilesCount + 1);

  std::vector<double> out;
  out.reserve(percentiles.size());
  for (uint64_t discrete_value : percentiles) {
    out.emplace_back(static_cast<double>(discrete_value) / kDiscreteMultiplier);
  }

  return out;
}

static std::vector<double> ParseSingle(const DataVector& path_counts,
                                       const DataVector& path_data,
                                       double* total, double* count) {
  CHECK(path_counts.size() == path_data.size());
  nc::DiscreteDistribution<uint64_t> dist;
  for (size_t i = 0; i < path_counts.size(); ++i) {
    double value = path_data[i].second;
    uint64_t value_discrete =
        static_cast<uint64_t>(kDiscreteMultiplier * value);

    double path_count;
    CHECK(std::modf(path_counts[i].second, &path_count) == 0.0);
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
  nc::viz::DataSeries2D aggregate_path_count;
};

static std::unique_ptr<Result> HandleSingleOptimizer(
    const std::string& optimizer) {
  std::string fields_regex = nc::StrCat(".*", optimizer, "$");

  LOG(INFO) << "Parsing metric " << kPathCountMetric << " fields "
            << fields_regex;
  DataMap path_counts = SimpleParseNumericData(
      FLAGS_metric_file, kPathCountMetric, fields_regex, 0, kMaxUint, 0);
  LOG(INFO) << "Done";

  auto result_ptr = nc::make_unique<Result>();
  LOG(INFO) << "Parsing metric " << kPathStretchRelMetric << " fields "
            << fields_regex;
  DataMap per_path_data = SimpleParseNumericData(
      FLAGS_metric_file, kPathStretchRelMetric, fields_regex, 0, kMaxUint, 0);
  DataMap per_path_penalty = SimpleParseNumericData(
      FLAGS_metric_file, kPathPenaltyMetric, fields_regex, 0, kMaxUint, 0);

  // Will plot the percentiles to get a CDF.
  std::vector<double> percentiles =
      Parse(path_counts, per_path_data, per_path_penalty, 100.0);
  for (size_t i = 0; i < percentiles.size(); ++i) {
    double p = percentiles[i];
    result_ptr->path_stretch_rel_data.data.emplace_back(p, i);
  }
  result_ptr->path_stretch_rel_data.label = optimizer;

  std::vector<double> per_aggregate_path_counts;
  DataMap per_aggregate_data =
      SimpleParseNumericData(FLAGS_metric_file, kAggregatePathCountMetric,
                             fields_regex, 0, kMaxUint, 0);
  for (const auto& key_and_data : per_aggregate_data) {
    const DataVector& path_data = key_and_data.second;
    for (const auto& datapoint : path_data) {
      per_aggregate_path_counts.emplace_back(datapoint.second);
    }
  }
  percentiles = nc::Percentiles(&per_aggregate_path_counts, kPercentilesCount);
  for (size_t i = 0; i < percentiles.size(); ++i) {
    double p = percentiles[i];
    result_ptr->aggregate_path_count.data.emplace_back(p, i);
  }

  return result_ptr;
}

static nc::viz::DataSeries1D FlattenData(const DataMap& data_map,
                                         const std::string& label) {
  std::vector<double> values;
  for (const auto& key_and_value : data_map) {
    for (const auto& v : key_and_value.second) {
      values.emplace_back(v.second);
    }
  }

  nc::viz::DataSeries1D data_series;
  data_series.data = std::move(values);
  data_series.label = label;

  return data_series;
}

static std::unique_ptr<nc::viz::DataSeries2D> HandleFocusOn(
    const std::string& optimizer) {
  std::string field_regex =
      nc::StrCat(".*", FLAGS_focus_on, ".*", optimizer, "$");
  LOG(INFO) << "Parsing " << kPathCountMetric << " for " << field_regex;
  DataMap path_count_data = SimpleParseNumericData(
      FLAGS_metric_file, kPathCountMetric, field_regex, 0, kMaxUint, 0);

  LOG(INFO) << "Parsing " << kPathStretchAbsMetric << " for " << field_regex;
  DataMap path_stretch_abs_data = SimpleParseNumericData(
      FLAGS_metric_file, kPathStretchAbsMetric, field_regex, 0, kMaxUint, 0);

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

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_metric_file.empty()) << "need --metric_file";

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

  std::vector<std::unique_ptr<Result>> result_ptrs =
      nc::RunInParallelWithResult<std::string, Result>(
          optimizers, [](const std::string& opt) {
            return HandleSingleOptimizer(opt);
          }, optimizers.size());

  PLOT("path_rel_stretch_out", result_ptrs[i]->path_stretch_rel_data);
  PLOT("aggregate_path_count_out", result_ptrs[i]->aggregate_path_count);

  DataMap runtime_data = SimpleParseNumericData(
      FLAGS_metric_file, "ctr_runtime_ms", ".*", 0, kMaxUint, 0);
  DataMap runtime_data_cached = SimpleParseNumericData(
      FLAGS_metric_file, "ctr_runtime_cached_ms", ".*", 0, kMaxUint, 0);
  nc::viz::PythonGrapher grapher("runtime_out");
  grapher.PlotCDF({},
                  {FlattenData(runtime_data, "Runtime (ms)"),
                   FlattenData(runtime_data_cached, "Runtime cached (ms)")});
}
