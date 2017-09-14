#include <gflags/gflags.h>
#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ncode_common/src/logging.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/viz/grapher.h"
#include "metrics/metrics_parser.h"

DEFINE_string(metrics_dir, "", "The metrics directory.");
DEFINE_string(optimizers, "MinMaxLD,CTRNFC,B4,CTR,MinMaxK10",
              "Optimizers to plot.");
DEFINE_bool(
    new_links, false,
    "Plot only data for new links, if false will plot for existing ones.");

static constexpr char kDecreasingDelayAggregateMetric[] =
    "decreasing_delay_aggregate_fraction";
static constexpr char kIncreasingDelayAggregateMetric[] =
    "increasing_delay_aggregate_fraction";
static constexpr char kOverloadDeltaMetric[] = "overload_delta_fraction";
static constexpr char kTotalDelayDeltaMetric[] = "total_delay_delta_fraction";

using namespace nc::metrics::parser;

using DataVector = std::vector<double>;
using DataMap = std::map<std::pair<std::string, std::string>, DataVector>;

static std::vector<double> DataForOptimizer(const std::string& opt,
                                            const std::string& metric) {
  std::string regex = nc::StrCat(".*:", opt);
  nc::StrAppend(&regex, FLAGS_new_links ? ":1$" : ":0$");
  DataMap runtime_data =
      SimpleParseNumericDataNoTimestamps(FLAGS_metrics_dir, metric, regex);

  std::vector<double> out;
  for (const auto& metric_and_data : runtime_data) {
    const std::vector<double>& data = metric_and_data.second;
    out.insert(out.end(), data.begin(), data.end());
  }

  return out;
}

static void PlotMetric(const std::string& metric) {
  std::vector<nc::viz::DataSeries1D> to_plot;
  std::vector<std::string> opts = nc::Split(FLAGS_optimizers, ",");
  for (const auto& opt : opts) {
    std::vector<double> data = DataForOptimizer(opt, metric);

    nc::viz::DataSeries1D data_series;
    data_series.data = std::move(data);
    data_series.label = opt;
    to_plot.emplace_back(data_series);
  }

  nc::viz::PythonGrapher grapher(nc::StrCat("top_grow_", metric));
  grapher.PlotCDF({}, to_plot);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_metrics_dir.empty()) << "need --metrics_dir";

  PlotMetric(kIncreasingDelayAggregateMetric);
  PlotMetric(kDecreasingDelayAggregateMetric);
  PlotMetric(kTotalDelayDeltaMetric);
  PlotMetric(kOverloadDeltaMetric);
}
