#include <gflags/gflags.h>
#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ncode_common/src/logging.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/substitute.h"
#include "ncode_common/src/viz/grapher.h"
#include "metrics/metrics_parser.h"

DEFINE_string(metrics_dir, "", "The metrics directory.");
DEFINE_string(optimizers, "MinMaxLD,CTRNFC,B4,CTR,MinMaxK10",
              "Optimizers to plot.");
DEFINE_string(change_threshold, "0\\.0",
              "Change threshold. Will match the string from the change "
              "threshold field in the metric.");
DEFINE_string(rank_links_for, "", "A topology/tm to rank links for");

static constexpr char kDecreasingDelayAggregateMetric[] =
    "decreasing_delay_aggregate_fraction";
static constexpr char kIncreasingDelayAggregateMetric[] =
    "increasing_delay_aggregate_fraction";
static constexpr char kDecreasingDelayFlowMetric[] =
    "decreasing_delay_flows_fraction";
static constexpr char kIncreasingDelayFlowMetric[] =
    "increasing_delay_flows_fraction";
static constexpr char kOverloadDeltaMetric[] = "overload_delta_fraction";
static constexpr char kTotalDelayDeltaMetric[] = "total_delay_delta_fraction";
static constexpr char kAbsoluteFlowStretchMetric[] = "absolute_path_stretch_ms";
static constexpr char kRelativeFlowStretchMetric[] = "relative_path_stretch";
static constexpr char kLinkDelayMetric[] = "link_delay_micros";

using namespace nc::metrics::parser;

using DataVector = std::vector<double>;
using DataMap = std::map<std::pair<std::string, std::string>, DataVector>;

using DistVector = std::vector<nc::DiscreteDistribution<int64_t>>;
using DistMap = std::map<std::pair<std::string, std::string>, DistVector>;

static std::vector<double> DataForOptimizer(const std::string& opt,
                                            const std::string& metric,
                                            bool new_links) {
  std::string regex = nc::Substitute(
      ".*:$0:$1(:$2)?", opt, std::to_string(new_links), FLAGS_change_threshold);
  nc::StrAppend(&regex, "$");

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
  std::string regex = nc::Substitute(
      ".*:$0:$1(:$2)?", opt, std::to_string(new_links), FLAGS_change_threshold);
  nc::StrAppend(&regex, "$");

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

  nc::viz::PythonGrapher grapher(nc::StrCat("top_grow_", metric, "_new_links_",
                                            new_links, "_thres_",
                                            FLAGS_change_threshold));
  nc::viz::PlotParameters1D params;
  params.title =
      nc::StrCat(title, " delay change at least ", FLAGS_change_threshold);
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
                  nc::viz::DataSeries2D, nc::viz::DataSeries2D>
PlotRankedLinks(const DataMap& increase_data, const DataMap& decrease_data,
                const DataMap& total_delta, const DataMap& link_delay_data,
                const std::vector<size_t>& order) {
  CHECK(increase_data.size() == 1);
  CHECK(decrease_data.size() == 1);
  CHECK(total_delta.size() == 1);
  CHECK(link_delay_data.size() == 1);

  const std::vector<double>& increase_v = increase_data.begin()->second;
  const std::vector<double>& decrease_v = decrease_data.begin()->second;
  const std::vector<double>& delta_v = total_delta.begin()->second;
  const std::vector<double>& link_delay_v = link_delay_data.begin()->second;
  CHECK(increase_v.size() == decrease_v.size());
  CHECK(delta_v.size() == decrease_v.size());
  CHECK(delta_v.size() == order.size());
  CHECK(link_delay_v.size() == order.size());

  nc::viz::DataSeries2D increase_to_plot;
  nc::viz::DataSeries2D decrease_to_plot;
  nc::viz::DataSeries2D total_to_plot;
  nc::viz::DataSeries2D link_delay_to_plot;
  for (size_t i = 0; i < order.size(); ++i) {
    size_t index = order[i];
    total_to_plot.data.emplace_back(i, delta_v[index]);
    increase_to_plot.data.emplace_back(i, increase_v[index]);
    decrease_to_plot.data.emplace_back(i, decrease_v[index]);
    link_delay_to_plot.data.emplace_back(i, link_delay_v[index]);
  }

  return std::make_tuple(increase_to_plot, decrease_to_plot, total_to_plot,
                         link_delay_to_plot);
}

std::vector<size_t> PickOrder(const std::string& field, bool new_links) {
  std::string regex =
      nc::Substitute("$0:$1:$2(:$3)?", field, "CTR", std::to_string(new_links),
                     FLAGS_change_threshold);
  nc::StrAppend(&regex, "$");
  DataMap total_ctr = SimpleParseNumericDataNoTimestamps(
      FLAGS_metrics_dir, kTotalDelayDeltaMetric, regex);

  regex = nc::Substitute("$0:$1:$2(:$3)?", field, "MinMaxLD",
                         std::to_string(new_links), FLAGS_change_threshold);
  DataMap total_minmax = SimpleParseNumericDataNoTimestamps(
      FLAGS_metrics_dir, kTotalDelayDeltaMetric, regex);

  CHECK(total_ctr.size() == 1);
  CHECK(total_minmax.size() == 1);

  const std::vector<double>& ctr_v = total_ctr.begin()->second;
  const std::vector<double>& minmax_v = total_minmax.begin()->second;
  CHECK(ctr_v.size() == minmax_v.size());

  std::vector<std::tuple<double, double, size_t>> all_data(ctr_v.size());
  for (size_t i = 0; i < ctr_v.size(); ++i) {
    all_data[i] = std::make_tuple(ctr_v[i], minmax_v[i], i);
  }
  std::sort(all_data.begin(), all_data.end());

  std::vector<size_t> indices(all_data.size());
  for (size_t i = 0; i < ctr_v.size(); ++i) {
    indices[i] = std::get<2>(all_data[i]);
  }

  return indices;
}

static void PlotLinks(const std::string& field, bool new_links) {
  std::vector<size_t> order = PickOrder(field, new_links);
  std::vector<std::string> opts = nc::Split(FLAGS_optimizers, ",");
  std::vector<nc::viz::DataSeries2D> increase_series;
  std::vector<nc::viz::DataSeries2D> decrease_series;
  std::vector<nc::viz::DataSeries2D> total_series;
  std::vector<nc::viz::DataSeries2D> link_delay_series;

  for (const auto& opt : opts) {
    std::string regex =
        nc::Substitute("$0:$1:$2(:$3)?", field, opt, std::to_string(new_links),
                       FLAGS_change_threshold);
    nc::StrAppend(&regex, "$");

    DataMap increase_data_map = SimpleParseNumericDataNoTimestamps(
        FLAGS_metrics_dir, kIncreasingDelayFlowMetric, regex);
    DataMap decrease_data_map = SimpleParseNumericDataNoTimestamps(
        FLAGS_metrics_dir, kDecreasingDelayFlowMetric, regex);
    DataMap delta_data_map = SimpleParseNumericDataNoTimestamps(
        FLAGS_metrics_dir, kTotalDelayDeltaMetric, regex);
    DataMap link_delay_data_map = SimpleParseNumericDataNoTimestamps(
        FLAGS_metrics_dir, kLinkDelayMetric, regex);

    nc::viz::DataSeries2D increase_data;
    nc::viz::DataSeries2D decrease_data;
    nc::viz::DataSeries2D total_data;
    nc::viz::DataSeries2D link_delay_data;
    std::tie(increase_data, decrease_data, total_data, link_delay_data) =
        PlotRankedLinks(increase_data_map, decrease_data_map, delta_data_map,
                        link_delay_data_map, order);

    increase_data.label = opt;
    decrease_data.label = opt;
    total_data.label = opt;
    link_delay_data.label = opt;

    increase_series.emplace_back(increase_data);
    decrease_series.emplace_back(decrease_data);
    total_series.emplace_back(total_data);
    link_delay_series.emplace_back(link_delay_data);
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

  nc::viz::PythonGrapher ld_grapher(
      nc::StrCat("top_grow_links_ranked_link_delay_new_links_", new_links));
  ld_grapher.PlotLine({}, link_delay_series);
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
