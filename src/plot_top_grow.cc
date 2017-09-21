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
DEFINE_string(rank_links_for, "", "A topology/tm to rank links for");

static constexpr char kDecreasingDelayAggregateMetric[] =
    "decreasing_delay_aggregate_fraction";
static constexpr char kIncreasingDelayAggregateMetric[] =
    "increasing_delay_aggregate_fraction";
static constexpr char kDecreasingDelayFlowMetric[] =
    "decreasing_delay_flow_fraction";
static constexpr char kIncreasingDelayFlowMetric[] =
    "increasing_delay_flow_fraction";
static constexpr char kOverloadDeltaMetric[] = "overload_delta_fraction";
static constexpr char kTotalDelayDeltaMetric[] = "total_delay_delta_fraction";
static constexpr char kAbsoluteFlowStretchMetric[] = "absolute_path_stretch_ms";
static constexpr char kRelativeFlowStretchMetric[] = "relative_path_stretch";

using namespace nc::metrics::parser;

using DataVector = std::vector<double>;
using DataMap = std::map<std::pair<std::string, std::string>, DataVector>;

using DistVector = std::vector<nc::DiscreteDistribution<int64_t>>;
using DistMap = std::map<std::pair<std::string, std::string>, DistVector>;

static std::vector<double> DataForOptimizer(const std::string& opt,
                                            const std::string& metric,
                                            bool new_links) {
  std::string regex = nc::StrCat(".*:", opt);
  nc::StrAppend(&regex, new_links ? ":1$" : ":0$");
  DataMap data_map =
      SimpleParseNumericDataNoTimestamps(FLAGS_metrics_dir, metric, regex);

  std::vector<double> out;
  for (const auto& metric_and_data : data_map) {
    const std::vector<double>& data = metric_and_data.second;
    out.insert(out.end(), data.begin(), data.end());
  }

  return out;
}

static nc::DiscreteDistribution<int64_t> DistDataForOptimizer(
    const std::string& opt, const std::string& metric, bool new_links) {
  std::string regex = nc::StrCat(".*:", opt);
  nc::StrAppend(&regex, new_links ? ":1$" : ":0$");
  DistMap data_map =
      SimpleParseDistributionDataNoTimestamps(FLAGS_metrics_dir, metric, regex);

  nc::DiscreteDistribution<int64_t> out;
  for (const auto& metric_and_data : data_map) {
    for (const nc::DiscreteDistribution<int64_t>& dist :
         metric_and_data.second) {
      out.Add(dist);
    }
  }

  return out;
}

static void PlotMetric(const std::string& metric, const std::string& title,
                       bool new_links) {
  std::vector<nc::viz::DataSeries1D> to_plot;
  std::vector<std::string> opts = nc::Split(FLAGS_optimizers, ",");
  for (const auto& opt : opts) {
    std::vector<double> data = DataForOptimizer(opt, metric, new_links);

    nc::viz::DataSeries1D data_series;
    data_series.data = std::move(data);
    data_series.label = opt;
    to_plot.emplace_back(data_series);
  }

  nc::viz::PythonGrapher grapher(
      nc::StrCat("top_grow_", metric, "_new_links_", new_links));
  nc::viz::PlotParameters1D params;
  params.title = title;
  grapher.PlotCDF(params, to_plot);
}

static void PlotStretch(const std::string& metric, double multiplier,
                        const std::string& title, bool new_links) {
  std::vector<nc::viz::DataSeries2D> to_plot;
  std::vector<std::string> opts = nc::Split(FLAGS_optimizers, ",");
  for (const auto& opt : opts) {
    nc::DiscreteDistribution<int64_t> dist =
        DistDataForOptimizer(opt, metric, new_links);

    std::vector<int64_t> percentiles = dist.Percentiles(10000);
    nc::viz::DataSeries2D data_series;
    for (size_t i = 0; i < percentiles.size(); ++i) {
      data_series.data.emplace_back(percentiles[i] * multiplier, i);
    }
    data_series.label = opt;
    to_plot.emplace_back(data_series);
  }

  nc::viz::PythonGrapher grapher(
      nc::StrCat("top_grow_", metric, "_new_links_", new_links));
  nc::viz::PlotParameters2D params;
  params.title = title;
  grapher.PlotLine(params, to_plot);
}

static std::tuple<nc::viz::DataSeries2D, nc::viz::DataSeries2D,
                  nc::viz::DataSeries2D>
PlotRankedLinks(const DataMap& increase_data, const DataMap& decrease_data,
                const DataMap& total_delta) {
  CHECK(increase_data.size() == 1);
  CHECK(decrease_data.size() == 1);
  CHECK(total_delta.size() == 1);

  const std::vector<double>& increase_v = increase_data.begin()->second;
  const std::vector<double>& decrease_v = decrease_data.begin()->second;
  const std::vector<double>& delta_v = total_delta.begin()->second;
  CHECK(increase_v.size() == decrease_v.size());
  CHECK(delta_v.size() == decrease_v.size());

  // Will rank all based on total_delta.
  std::vector<std::tuple<double, double, double>> all_data;
  for (size_t i = 0; i < delta_v.size(); ++i) {
    all_data.emplace_back(delta_v[i], increase_v[i], decrease_v[i]);
  }
  std::sort(all_data.begin(), all_data.end());

  nc::viz::DataSeries2D increase_to_plot;
  nc::viz::DataSeries2D decrease_to_plot;
  nc::viz::DataSeries2D total_to_plot;
  for (size_t i = 0; i < delta_v.size(); ++i) {
    total_to_plot.data.emplace_back(i, std::get<0>(all_data[i]));
    increase_to_plot.data.emplace_back(i, std::get<1>(all_data[i]));
    decrease_to_plot.data.emplace_back(i, std::get<2>(all_data[i]));
  }

  return std::make_tuple(increase_to_plot, decrease_to_plot, total_to_plot);
}

static void PlotLinks(const std::string& field, bool new_links) {
  std::vector<std::string> opts = nc::Split(FLAGS_optimizers, ",");
  std::vector<nc::viz::DataSeries2D> increase_series;
  std::vector<nc::viz::DataSeries2D> decrease_series;
  std::vector<nc::viz::DataSeries2D> total_series;

  for (const auto& opt : opts) {
    std::string regex = nc::StrCat(field, ".*:", opt);
    nc::StrAppend(&regex, new_links ? ":1$" : ":0$");
    DataMap increase_data_map = SimpleParseNumericDataNoTimestamps(
        FLAGS_metrics_dir, kIncreasingDelayFlowMetric, regex);
    DataMap decrease_data_map = SimpleParseNumericDataNoTimestamps(
        FLAGS_metrics_dir, kDecreasingDelayFlowMetric, regex);
    DataMap delta_data_map = SimpleParseNumericDataNoTimestamps(
        FLAGS_metrics_dir, kTotalDelayDeltaMetric, regex);

    nc::viz::DataSeries2D increase_data;
    nc::viz::DataSeries2D decrease_data;
    nc::viz::DataSeries2D total_data;
    std::tie(increase_data, decrease_data, total_data) =
        PlotRankedLinks(increase_data_map, decrease_data_map, delta_data_map);

    increase_data.label = opt;
    decrease_data.label = opt;
    total_data.label = opt;

    increase_series.emplace_back(increase_data);
    decrease_series.emplace_back(decrease_data);
    total_series.emplace_back(total_data);
  }

  nc::viz::PythonGrapher increase_grapher(
      nc::StrCat("top_grow_links_ranked_increase_new_links_", new_links));
  increase_grapher.PlotLine({}, increase_series);

  nc::viz::PythonGrapher decrease_grapher(
      nc::StrCat("top_grow_links_ranked_decrease_new_links_", new_links));
  decrease_grapher.PlotLine({}, decrease_series);

  nc::viz::PythonGrapher total_grapher(
      nc::StrCat("top_grow_links_ranked_total_new_links_", new_links));
  total_grapher.PlotLine({}, total_series);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_metrics_dir.empty()) << "need --metrics_dir";

  if (!FLAGS_rank_links_for.empty()) {
    PlotLinks(FLAGS_rank_links_for, true);
    PlotLinks(FLAGS_rank_links_for, false);
    return 0;
  }

  PlotMetric(kIncreasingDelayAggregateMetric,
             "Fraction of aggregates with increasing delay (new links)", true);
  PlotMetric(kDecreasingDelayAggregateMetric,
             "Fraction of aggregates with decreasing delay (new links)", true);
  PlotMetric(kIncreasingDelayFlowMetric,
             "Fraction of flows with increasing delay (new links)", true);
  PlotMetric(kDecreasingDelayFlowMetric,
             "Fraction of flows with decreasing delay (new links)", true);
  PlotMetric(kTotalDelayDeltaMetric, "Total delay delta (new links)", true);
  PlotMetric(kOverloadDeltaMetric, "Overload (new links)", true);
  PlotStretch(kAbsoluteFlowStretchMetric, 1.0,
              "Absolute flow stretch (new links)", true);
  PlotStretch(kRelativeFlowStretchMetric, 1 / 10000.0,
              "Relative flow stretch (new links)", true);

  PlotMetric(kIncreasingDelayAggregateMetric,
             "Fraction of aggregates with increasing delay (old links)", false);
  PlotMetric(kDecreasingDelayAggregateMetric,
             "Fraction of aggregates with decreasing delay (old links)", false);
  PlotMetric(kIncreasingDelayFlowMetric,
             "Fraction of flows with increasing delay (old links)", false);
  PlotMetric(kDecreasingDelayFlowMetric,
             "Fraction of flows with decreasing delay (old links)", false);
  PlotMetric(kTotalDelayDeltaMetric, "Total delay delta (old links)", false);
  PlotMetric(kOverloadDeltaMetric, "Overload (old links)", false);
  PlotStretch(kAbsoluteFlowStretchMetric, 1.0,
              "Absolute flow stretch (old links)", false);
  PlotStretch(kRelativeFlowStretchMetric, 1 / 10000.0,
              "Relative flow stretch (old links)", false);
}
