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

DEFINE_string(sp_input, "", "The sp metrics file.");
DEFINE_string(minmax_input, "", "The minmax metrics file.");
DEFINE_string(metric, "tcp_source_completion_times",
              "The completion times metric.");
DEFINE_string(count_metric, "tcp_source_close_count",
              "The metric that identifies objects.");

using namespace nc;
using namespace nc::metrics::parser;

using DataVector = std::vector<std::pair<uint64_t, double>>;

static std::vector<std::chrono::milliseconds> GetDeltas(
    const DataVector& sp_data, const DataVector& other_data,
    const DataVector& sp_count_data, const DataVector& other_count_data) {
  size_t min_data_size = std::min(sp_data.size(), other_data.size());
  std::vector<std::chrono::milliseconds> out;
  for (size_t i = 0; i < min_data_size; ++i) {
    const auto& d1 = sp_data[i];
    const auto& d2 = other_data[i];

    if (other_data[i].second < sp_data[i].second) {
      LOG(ERROR) << other_data[i].first << " vs " << sp_data[i].first << " "
                 << (i + 1) << "/" << min_data_size << " SP completion time "
                 << sp_count_data[i].second << " other completion time "
                 << other_count_data[i].second;
      continue;
    }
    uint64_t delta = other_data[i].second - sp_data[i].second;
    out.emplace_back(std::chrono::milliseconds(delta));
  }

  return out;
}

static void PlotCompletionTimes(
    const std::vector<std::chrono::milliseconds>& minmax_data) {
  std::vector<double> data_ms;
  for (auto v : minmax_data) {
    data_ms.emplace_back(v.count());
  }

  nc::viz::DataSeries1D data_series;
  data_series.data = std::move(data_ms);
  data_series.label = "MinMax";

  nc::viz::PythonGrapher grapher("ct_plot_out");
  grapher.PlotCDF({}, {data_series});
  LOG(INFO) << "Plotted " << data_ms.size() << " values";
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_sp_input.empty()) << "need --sp_input";
  CHECK(!FLAGS_minmax_input.empty()) << "need --sp_input";

  std::map<std::pair<std::string, std::string>, DataVector> sp_data =
      SimpleParseNumericData(FLAGS_sp_input, FLAGS_metric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);

  std::map<std::pair<std::string, std::string>, DataVector> minmax_data =
      SimpleParseNumericData(FLAGS_minmax_input, FLAGS_metric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);
  CHECK(sp_data.size() == minmax_data.size());

  std::map<std::pair<std::string, std::string>, DataVector> sp_count_data =
      SimpleParseNumericData(FLAGS_sp_input, "tcp_source_close_count", ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);

  std::map<std::pair<std::string, std::string>, DataVector> minmax_count_data =
      SimpleParseNumericData(FLAGS_minmax_input, "tcp_source_close_count", ".*",
                             0, std::numeric_limits<uint64_t>::max(), 0);

  std::vector<std::chrono::milliseconds> all_values;
  for (auto& id_and_data : sp_data) {
    const std::pair<std::string, std::string>& id = id_and_data.first;
    DataVector& data = id_and_data.second;
    DataVector& other_data = nc::FindOrDieNoPrint(minmax_data, id);
    DataVector& count_data =
        nc::FindOrDieNoPrint(sp_count_data, {FLAGS_count_metric, id.second});
    DataVector& other_count_data = nc::FindOrDieNoPrint(
        minmax_count_data, {FLAGS_count_metric, id.second});
    CHECK(id.first == FLAGS_metric);

    std::vector<std::chrono::milliseconds> values =
        GetDeltas(data, other_data, count_data, other_count_data);
    all_values.insert(all_values.end(), values.begin(), values.end());
  }

  PlotCompletionTimes(all_values);
}
