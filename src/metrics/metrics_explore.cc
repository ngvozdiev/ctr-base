#include <gflags/gflags.h>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "ncode/common.h"
#include "ncode/logging.h"
#include "ncode/stats.h"
#include "ncode/strutil.h"
#include "ncode/substitute.h"
#include "ncode/viz/grapher.h"
#include "metrics_parser.h"

DEFINE_string(input, "", "The metrics file.");
DEFINE_string(metric, "",
              "The id of the metric. If missing will print summary of all "
              "metrics. Can be a regex and can match more than one metric.");
DEFINE_string(fields, ".*", "The fields to match.");
DEFINE_string(output, "metrics_explorer_out", "Output directory for plots");
DEFINE_string(plot_type, "", "Either cdf or line, no plot generated if empty.");
DEFINE_bool(plot_combined, false,
            "If true will combine vales from all fields sets and only "
            "produce one plot.");
DEFINE_bool(plot_limiting, false,
            "If true will only plot the last value of each field set");
DEFINE_bool(incremental, false,
            "If true the metric(s) matched will be assumed to be incrementally "
            "increasing.");
DEFINE_uint64(line_plot_bin_size, 1,
              "Will bin every N values based on "
              "timestamp and plot the means in each bin.");
DEFINE_uint64(plot_timestamp_min, 0, "Min timestamp to plot");
DEFINE_uint64(plot_timestamp_max, std::numeric_limits<uint64_t>::max(),
              "Max timestamp to plot");

using namespace nc;
using namespace nc::metrics::parser;

// For a series (t1, v1), (t2, v2), (t3, v3) ... returns (t2, v2 - v1), (t3, v3
// - v2), ...
static std::vector<std::pair<uint64_t, double>> TakeDeltas(
    const std::vector<std::pair<uint64_t, double>>& values) {
  CHECK(!values.empty());

  std::vector<std::pair<uint64_t, double>> out;
  for (size_t i = 0; i < values.size() - 1; ++i) {
    uint64_t timestamp = values[i + 1].first;
    double value = values[i + 1].second - values[i].second;
    out.emplace_back(timestamp, value);
  }

  return out;
}

// Plots a CDF of the values in one (or more) metrics. The input map is as
// returned by any of the query functions in metrics_parser.
static void PlotCDF(
    const std::map<std::pair<std::string, std::string>,
                   std::vector<std::pair<uint64_t, double>>>& data) {
  if (data.empty()) {
    LOG(ERROR) << "No numeric data";
    return;
  }

  std::vector<viz::DataSeries1D> all_data_to_plot;
  if (FLAGS_plot_combined) {
    viz::DataSeries1D to_plot;
    std::vector<double> all_values;
    for (const auto& label_and_values : data) {
      const std::vector<std::pair<uint64_t, double>>& values =
          label_and_values.second;
      for (const auto& timestamp_and_value : values) {
        double value = timestamp_and_value.second;
        all_values.emplace_back(value);
      }
    }

    to_plot.data = std::move(all_values);
    all_data_to_plot.emplace_back(std::move(to_plot));
  } else {
    for (const auto& label_and_values : data) {
      const std::string& label_id = label_and_values.first.second;
      const std::vector<std::pair<uint64_t, double>>& values =
          label_and_values.second;

      viz::DataSeries1D to_plot;
      for (const auto& timestamp_and_value : values) {
        double value = timestamp_and_value.second;
        to_plot.data.emplace_back(value);
      }

      to_plot.label = label_id;
      all_data_to_plot.emplace_back(std::move(to_plot));
    }
  }

  viz::CDFPlot plot;
  plot.AddData(all_data_to_plot);
  plot.PlotToDir(FLAGS_output);
  LOG(INFO) << "Saved script to plot data at " << FLAGS_output;
}

static void PlotTimeSeries(
    const std::map<std::pair<std::string, std::string>,
                   std::vector<std::pair<uint64_t, double>>>& data) {
  if (data.empty()) {
    LOG(ERROR) << "No numeric data";
    return;
  }

  std::vector<viz::DataSeries2D> all_data_to_plot;
  for (const auto& label_and_values : data) {
    const std::string& label_id = label_and_values.first.second;
    const std::vector<std::pair<uint64_t, double>>& values =
        label_and_values.second;

    viz::DataSeries2D to_plot;
    to_plot.label = label_id;

    for (const auto& timestamp_and_value : values) {
      uint64_t timestamp = timestamp_and_value.first;
      double value = timestamp_and_value.second;
      to_plot.data.emplace_back(timestamp, value);
    }
    all_data_to_plot.emplace_back(std::move(to_plot));
  }

  viz::PlotParameters2D plot_params;
  plot_params.x_bin_size = FLAGS_line_plot_bin_size;
  viz::LinePlot plot(plot_params);
  plot.AddData(all_data_to_plot);
  plot.PlotToDir(FLAGS_output);

  LOG(INFO) << "Saved script to plot data at " << FLAGS_output;
}

static void HandlePlot(const std::string& input_file) {
  std::map<std::pair<std::string, std::string>,
           std::vector<std::pair<uint64_t, double>>> data =
      SimpleParseNumericData(
          input_file, FLAGS_metric, FLAGS_fields, FLAGS_plot_timestamp_min,
          FLAGS_plot_timestamp_max,
          FLAGS_plot_limiting ? std::numeric_limits<uint64_t>::max() : 0);

  if (FLAGS_incremental) {
    for (auto& id_and_data : data) {
      std::vector<std::pair<uint64_t, double>>& metric_data =
          id_and_data.second;
      metric_data = TakeDeltas(metric_data);
    }
  }

  if (FLAGS_plot_type == "cdf") {
    PlotCDF(data);
  } else if (FLAGS_plot_type == "line") {
    PlotTimeSeries(data);
  }
}

static void ProcessSingleFile(const std::string& input_file) {
  if (!FLAGS_metric.empty()) {
    if (!FLAGS_plot_type.empty()) {
      HandlePlot(input_file);
      return;
    }

    std::map<std::pair<std::string, std::string>,
             std::vector<std::pair<uint64_t, double>>> data =
        SimpleParseNumericData(input_file, FLAGS_metric, FLAGS_fields, 0,
                               std::numeric_limits<uint64_t>::max(), 0);
    std::cout << DataSummaryToString(data);
    return;
  }

  MetricsParser parser(input_file);
  Manifest manifest = parser.ParseManifest();
  std::cout << manifest.FullToString();
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_input.empty()) << "need --input";

  std::vector<std::string> inputs = Glob(FLAGS_input);
  for (const auto& input : inputs) {
    ProcessSingleFile(input);
  }
}
