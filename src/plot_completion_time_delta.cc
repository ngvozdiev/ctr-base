#include <gflags/gflags.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/viz/grapher.h"
#include "metrics/metrics_parser.h"

DEFINE_string(sp_input, "", "The shortest path metrics file.");
DEFINE_string(other_inputs, "",
              "A comma-separated list of other metric files to examine.");
DEFINE_string(labels, "",
              "A comma-separated list of labels. "
              "One per metric file in other_inputs");
DEFINE_string(metric, "tcp_source_completion_times",
              "The completion times metric.");
DEFINE_bool(absolute, false,
            "If true will plot CDFs of absolute completion times, if false "
            "will plot CDFs of CTs relative to those of the shortest path.");
DEFINE_uint64(grace_period_sec, 60,
              "Will ignore completion times that are generated within the "
              "first this many seconds.");

using namespace nc;
using namespace nc::metrics::parser;

using DataVector = std::vector<std::pair<uint64_t, double>>;

// Will assume that the timestamps are in picoseconds.
static constexpr uint64_t kTimestampMultiplier = 1000000000000UL;

static std::vector<std::chrono::milliseconds> GetDeltas(
    const DataVector& sp_data, const DataVector& other_data) {
  uint64_t threshold_timestamp = kTimestampMultiplier * FLAGS_grace_period_sec;

  size_t min_data_size = std::min(sp_data.size(), other_data.size());
  std::vector<std::chrono::milliseconds> out;
  for (size_t i = 0; i < min_data_size; ++i) {
    uint64_t time = std::max(other_data[i].first, sp_data[i].first);
    if (time < threshold_timestamp) {
      continue;
    }

    if (other_data[i].second < sp_data[i].second) {
      LOG(ERROR) << other_data[i].first << " vs " << sp_data[i].first << " "
                 << (i + 1) << "/" << min_data_size << " SP completion time "
                 << sp_data[i].second << " other completion time "
                 << other_data[i].second;
      out.emplace_back(std::chrono::milliseconds(0));
      continue;
    }
    uint64_t delta = other_data[i].second - sp_data[i].second;
    out.emplace_back(std::chrono::milliseconds(delta));
  }

  return out;
}

static nc::viz::DataSeries1D GetDeltasFromSP(
    const std::map<std::pair<std::string, std::string>, DataVector>& sp_data,
    const std::string& file, const std::string& label) {
  nc::viz::DataSeries1D data_series;
  data_series.label = label;

  std::map<std::pair<std::string, std::string>, DataVector> other_data =
      SimpleParseNumericData(file, FLAGS_metric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);
  CHECK(sp_data.size() == other_data.size()) << sp_data.size() << " vs "
                                             << other_data.size();
  for (const auto& id_and_data : sp_data) {
    const std::pair<std::string, std::string>& id = id_and_data.first;
    const DataVector& data = id_and_data.second;
    const DataVector& other_data_inner = nc::FindOrDieNoPrint(other_data, id);
    CHECK(id.first == FLAGS_metric);

    std::vector<std::chrono::milliseconds> values =
        GetDeltas(data, other_data_inner);
    for (const auto& value : values) {
      data_series.data.emplace_back(value.count());
    }
  }

  return data_series;
}

static std::vector<std::string> GetInputFiles() {
  std::vector<std::string> split = nc::Split(FLAGS_other_inputs, ",");
  std::vector<std::string> out;
  for (const std::string& v : split) {
    std::vector<std::string> globbed = nc::Glob(v);
    out.insert(out.end(), globbed.begin(), globbed.end());
  }
  return out;
}

static std::vector<std::string> GetLabels() {
  return nc::Split(FLAGS_labels, ",");
}

nc::viz::DataSeries1D GetAbsoluteDataFromFile(const std::string& file,
                                              const std::string& label) {
  nc::viz::DataSeries1D data_series;
  data_series.label = label;

  std::map<std::pair<std::string, std::string>, DataVector> data =
      SimpleParseNumericData(file, FLAGS_metric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);
  for (const auto& id_and_data : data) {
    const DataVector& data_vector = id_and_data.second;
    for (const auto& value : data_vector) {
      data_series.data.emplace_back(value.second);
    }
  }
  return data_series;
}

static void HandleAbsolute(const std::vector<std::string>& input_files,
                           const std::vector<std::string>& labels) {
  std::vector<nc::viz::DataSeries1D> to_plot = {
      GetAbsoluteDataFromFile(FLAGS_sp_input, "SP")};
  for (size_t i = 0; i < input_files.size(); ++i) {
    nc::viz::DataSeries1D data_series =
        GetAbsoluteDataFromFile(input_files[i], labels[i]);
    to_plot.emplace_back(data_series);
  }

  nc::viz::CDFPlot plot;
  plot.AddData(to_plot);
  plot.PlotToDir("ct_plot_out");
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_sp_input.empty()) << "need --sp_input";

  std::map<std::pair<std::string, std::string>, DataVector> sp_data =
      SimpleParseNumericData(FLAGS_sp_input, FLAGS_metric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);

  std::vector<std::string> input_files = GetInputFiles();
  CHECK(!input_files.empty()) << "need at least one file in --other_inputs";
  std::vector<std::string> labels = GetLabels();
  CHECK(labels.size() == input_files.size());
  if (FLAGS_absolute) {
    HandleAbsolute(input_files, labels);
    return 0;
  }

  std::vector<nc::viz::DataSeries1D> to_plot;
  for (size_t i = 0; i < input_files.size(); ++i) {
    nc::viz::DataSeries1D data_series =
        GetDeltasFromSP(sp_data, input_files[i], labels[i]);
    to_plot.emplace_back(data_series);
  }

  nc::viz::CDFPlot plot;
  plot.AddData(to_plot);
  plot.PlotToDir("ct_plot_out");
  return 0;
}
