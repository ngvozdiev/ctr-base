#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/logging.h"
#include "ncode_common/src/stats.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/substitute.h"
#include "ncode_common/src/viz/grapher.h"
#include "metrics/metrics_parser.h"

DEFINE_string(input, "", "The metrics file.");
DEFINE_double(at_link_multiplier, 0.9, "What headroom value to plot.");
DEFINE_string(output_pattern, "headroom_delay_plot_out_$0",
              "What the output folder should be named. $0 will be replaced "
              "with 'at_link_multiplier'");

static constexpr char kStrideSizeMetric[] = "stride_size";
static constexpr char kTotalDelayMetric[] = "total_delay_fraction";

using namespace nc;
using namespace nc::metrics::parser;

using NumDataVector = std::vector<double>;
using StrPair = std::pair<std::string, std::string>;

// Returns the delay fraction at a given headroom (or the closest to the given
// headroom).
static double GetValueAtLinkMulitplier(const std::vector<double>& values,
                                       double stride) {
  // Link multiplier starts at 1.0, and decreases by 'stride' per item in
  // values. 1 - stride * i = multiplier; i = (1 - multiplier) / stride
  CHECK(FLAGS_at_link_multiplier < 1);
  size_t index = (1 - FLAGS_at_link_multiplier) / stride;
  if (index >= values.size()) {
    return std::numeric_limits<double>::max();
  }

  return values[index];
}

static std::vector<double> GetValuesAtLinkMultiplier(
    const std::vector<std::vector<double>>& values, double stride) {
  std::vector<double> out;
  for (const auto& value_vector : values) {
    double v = GetValueAtLinkMulitplier(value_vector, stride);
    if (v != std::numeric_limits<double>::max()) {
      out.emplace_back(v);
    }
  }

  return out;
}

static std::map<std::string, std::vector<std::vector<double>>> Group() {
  std::map<StrPair, NumDataVector> data =
      SimpleParseNumericDataNoTimestamps(FLAGS_input, kTotalDelayMetric, ".*");

  std::map<std::string, std::vector<std::vector<double>>> out;
  for (const auto& id_and_data : data) {
    const StrPair& id = id_and_data.first;
    const NumDataVector& data = id_and_data.second;

    std::vector<std::string> id_split = nc::Split(id.second, ":");
    CHECK(id_split.size() == 2);
    const std::string& topology = id_split[0];

    out[topology].emplace_back(data);
  }

  return out;
}

static double GetStride() {
  std::map<StrPair, NumDataVector> data =
      SimpleParseNumericDataNoTimestamps(FLAGS_input, kStrideSizeMetric, ".*");
  CHECK(data.size() == 1);
  CHECK(data.begin()->second.size() == 1);
  return data.begin()->second.front();
}

static std::map<std::string, std::vector<double>> GroupByLinkMultiplier() {
  double stride = GetStride();
  std::map<std::string, std::vector<std::vector<double>>> grouped = Group();
  std::map<std::string, std::vector<double>> out;
  for (const auto& topology_and_values : grouped) {
    const std::vector<std::vector<double>>& values = topology_and_values.second;
    out[topology_and_values.first] = GetValuesAtLinkMultiplier(values, stride);
  }

  return out;
}

static void PlotValues() {
  std::map<std::string, std::vector<double>> by_multiplier =
      GroupByLinkMultiplier();

  // Median, min, max for each topology.
  std::vector<std::tuple<double, double, double>> values;
  for (const auto& topology_and_data : by_multiplier) {
    std::vector<double> data = topology_and_data.second;
    if (data.empty()) {
      continue;
    }

    std::vector<double> p = nc::Percentiles(&data, 100);
    values.emplace_back(p[50], p[0], p[100]);
  }

  std::sort(values.begin(), values.end());
  std::vector<std::pair<double, double>> y_median;
  std::vector<std::pair<double, double>> y_min;
  std::vector<std::pair<double, double>> y_max;

  for (size_t i = 0; i < values.size(); ++i) {
    y_min.emplace_back(i, std::get<1>(values[i]));
    y_median.emplace_back(i, std::get<0>(values[i]));
    y_max.emplace_back(i, std::get<2>(values[i]));
  }

  nc::viz::PythonGrapher grapher(
      nc::Substitute(FLAGS_output_pattern.c_str(), FLAGS_at_link_multiplier));
  nc::viz::PlotParameters2D params;
  params.title = nc::StrCat("Total delay multiplier at link multiplier ",
                            FLAGS_at_link_multiplier);
  params.x_label = "topology (ranked by median delay multiplier)";
  params.y_label = "total delay multiplier";
  grapher.PlotLine(params, {{"median", y_median}, {"max", y_max}});
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_input.empty()) << "need --input";

  PlotValues();
  return 0;
}
