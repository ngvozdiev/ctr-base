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
DEFINE_string(minmax_input, "", "The MinMax metrics file.");
DEFINE_string(ctr_input, "", "The CTR metrics file.");
DEFINE_string(metric, "tcp_source_completion_times",
              "The completion times metric.");

using namespace nc;
using namespace nc::metrics::parser;

using DataVector = std::vector<std::pair<uint64_t, double>>;

static std::vector<std::chrono::milliseconds> GetDeltas(
    const DataVector& sp_data, const DataVector& other_data) {
  size_t min_data_size = std::min(sp_data.size(), other_data.size());
  std::vector<std::chrono::milliseconds> out;
  for (size_t i = 0; i < min_data_size; ++i) {
    if (other_data[i].second < sp_data[i].second) {
      LOG(ERROR) << other_data[i].first << " vs " << sp_data[i].first << " "
                 << (i + 1) << "/" << min_data_size << " SP completion time "
                 << sp_data[i].second << " other completion time "
                 << other_data[i].second;
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
  CHECK(sp_data.size() == other_data.size());
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

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_sp_input.empty()) << "need --sp_input";
  CHECK(!FLAGS_minmax_input.empty()) << "need --sp_input";

  std::map<std::pair<std::string, std::string>, DataVector> sp_data =
      SimpleParseNumericData(FLAGS_sp_input, FLAGS_metric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);

  nc::viz::DataSeries1D min_max_data_series =
      GetDeltasFromSP(sp_data, FLAGS_minmax_input, "MinMax");
  nc::viz::DataSeries1D ctr_data_series =
      GetDeltasFromSP(sp_data, FLAGS_ctr_input, "CTR");

  nc::viz::PythonGrapher grapher("ct_plot_out");
  grapher.PlotCDF({}, {min_max_data_series, ctr_data_series});
  LOG(INFO) << "Plotted " << ctr_data_series.data.size() << " values";
}
