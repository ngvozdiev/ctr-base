#include <gflags/gflags.h>
#include <ncode/logging.h>
#include <ncode/map_util.h>
#include <ncode/net/net_common.h>
#include <ncode/strutil.h>
#include <ncode/viz/grapher.h>
#include <stddef.h>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "metrics/metrics_parser.h"

DEFINE_string(input, "", "The metrics file.");
DEFINE_string(input_pinned, "", "The metrics file for a run with pinned means");
DEFINE_uint64(line_plot_bin_size, 1,
              "Will bin every N values based on "
              "timestamp and plot the means in each bin.");
DEFINE_uint64(n, 30, "Will plot the top N links");

static constexpr char kLinkUtilizationMetric[] = "link_utilization";
static constexpr char kLinkRateMetric[] = "link_rate_Mbps";
static constexpr char kBinSizeMetric[] = "bin_size_ms";
static constexpr char kQueueSizeMetric[] = "queue_size";
static constexpr char kRouteAddMetric[] = "route_add";
static constexpr char kRouteRemoveMetric[] = "route_remove";
static constexpr char kRouteUpdateMetric[] = "route_update";
static constexpr char kPropDelayMetric[] = "prop_delay_per_packet";

using namespace nc;
using namespace nc::metrics::parser;

using NumDataVector = std::vector<std::pair<uint64_t, double>>;
using StrPair = std::pair<std::string, std::string>;

static std::chrono::milliseconds GetBinSize() {
  std::map<StrPair, std::vector<double>> data =
      SimpleParseNumericDataNoTimestamps(FLAGS_input, kBinSizeMetric, ".*");
  CHECK(data.size() == 1);
  const std::vector<double>& values = data.begin()->second;
  CHECK(values.size() == 1);
  return std::chrono::milliseconds(static_cast<uint64_t>(values[0]));
}

static nc::net::Bandwidth GetLinkRate(const std::string& link_src_and_dst) {
  std::map<StrPair, std::vector<double>> data =
      SimpleParseNumericDataNoTimestamps(FLAGS_input, kLinkRateMetric, ".*");
  const std::vector<double>& values = nc::FindOrDieNoPrint(
      data, std::make_pair(kLinkRateMetric, link_src_and_dst));
  CHECK(values.size() == 1);
  return nc::net::Bandwidth::FromMBitsPerSecond(values[0]);
}

static void ParseSingleUpdateState(
    const std::string& metric, std::map<uint64_t, std::vector<uint64_t>>* map) {
  std::map<StrPair, NumDataVector> data = SimpleParseNumericData(
      FLAGS_input, metric, ".*", 0, std::numeric_limits<uint64_t>::max(), 0);
  CHECK(data.size() == 1);
  const NumDataVector& data_vector = data.begin()->second;
  for (const auto& timestamp_and_value : data_vector) {
    (*map)[timestamp_and_value.first].emplace_back(timestamp_and_value.second);
  }
}

static void PrintUpdateStats() {
  std::map<uint64_t, std::vector<uint64_t>> counts;
  ParseSingleUpdateState(kRouteAddMetric, &counts);
  ParseSingleUpdateState(kRouteRemoveMetric, &counts);
  ParseSingleUpdateState(kRouteUpdateMetric, &counts);

  for (const auto& timestamp_and_counts : counts) {
    uint64_t timestamp = timestamp_and_counts.first;
    const std::vector<uint64_t>& counts = timestamp_and_counts.second;
    CHECK(counts.size() == 3);
    LOG(INFO) << "Opt at " << timestamp << " add " << counts[0] << " remove "
              << counts[1] << " update " << counts[2];
  }
}

static double GetAverageUtilization(const std::string& link_src_and_dst,
                                    const NumDataVector& link_utilization,
                                    std::chrono::milliseconds bin_size) {
  nc::net::Bandwidth link_rate = GetLinkRate(link_src_and_dst);

  double total_bytes = 0;
  for (const auto& utilization : link_utilization) {
    total_bytes += utilization.second;
  }

  double bytes_per_second = link_rate.bps() / 8.0;
  double bins_per_second = 1000.0 / bin_size.count();
  double bytes_per_bin = bytes_per_second / bins_per_second;

  size_t bin_count = link_utilization.size();
  double total_link_bytes = bytes_per_bin * bin_count;
  return total_bytes / total_link_bytes;
}

static std::chrono::milliseconds GetMaxQueueSize(
    const std::string& link_src_and_dst,
    const NumDataVector& link_queue_sizes) {
  double max_queue_size_bytes = 0;
  for (const auto& queue_sizes : link_queue_sizes) {
    max_queue_size_bytes = std::max(queue_sizes.second, max_queue_size_bytes);
  }

  nc::net::Bandwidth link_rate = GetLinkRate(link_src_and_dst);
  double bytes_per_sec = link_rate.bps() / 8.0;
  double queue_size_sec = max_queue_size_bytes / bytes_per_sec;
  return std::chrono::milliseconds(
      static_cast<uint64_t>(queue_size_sec * 1000));
}

static void PlotLink(const std::string& link_src_and_dst,
                     const NumDataVector& link_utilization,
                     const std::string& output,
                     std::chrono::milliseconds bin_size) {
  nc::net::Bandwidth link_rate = GetLinkRate(link_src_and_dst);

  viz::PlotParameters2D plot_params(link_src_and_dst, "time", "rate (Mbps)");
  plot_params.x_bin_size = FLAGS_line_plot_bin_size;
  viz::LinePlot plot(plot_params);

  std::vector<double> y_values;
  std::vector<double> x_values;
  std::vector<double> link_rate_y_values;
  double bins_per_second = 1000.0 / bin_size.count();
  for (const auto& utilization : link_utilization) {
    double bits_in_bin = utilization.second * 8;
    double rate_bps = bits_in_bin * bins_per_second;
    y_values.emplace_back(rate_bps / 1000000.0);
    link_rate_y_values.emplace_back(link_rate.Mbps());
    x_values.emplace_back(utilization.first);
  }
  plot.AddData("traffic", x_values, y_values);
  plot.AddData("link rate", x_values, link_rate_y_values);
  plot.PlotToDir(output);
}

static void PlotTopNLinks(size_t n) {
  std::chrono::milliseconds bin_size = GetBinSize();

  std::map<StrPair, NumDataVector> data_pinned =
      SimpleParseNumericData(FLAGS_input_pinned, kLinkUtilizationMetric, ".*",
                             0, std::numeric_limits<uint64_t>::max(), 0);
  std::map<StrPair, NumDataVector> data =
      SimpleParseNumericData(FLAGS_input, kLinkUtilizationMetric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);

  std::map<double, std::map<StrPair, NumDataVector>::const_iterator> top_links;
  for (auto it = data_pinned.cbegin(); it != data_pinned.cend(); ++it) {
    const std::string& src_and_dst = it->first.second;
    const NumDataVector& data_vector = it->second;

    double avg_utilization =
        GetAverageUtilization(src_and_dst, data_vector, bin_size);
    top_links.emplace(avg_utilization, it);
  }

  uint32_t top_n = 0;
  for (auto rit = top_links.rbegin(); rit != top_links.rend(); ++rit) {
    double avg_utilization = rit->first;
    auto data_it = rit->second;
    const std::string& src_and_dst = data_it->first.second;
    LOG(INFO) << top_n << " " << src_and_dst << " " << avg_utilization;

    const NumDataVector& data_vector_pinned = data_it->second;
    const NumDataVector& data_vector_reg =
        nc::FindOrDieNoPrint(data, data_it->first);

    PlotLink(src_and_dst, data_vector_reg,
             nc::StrCat("link_top_", top_n, "_reg"), bin_size);
    PlotLink(src_and_dst, data_vector_pinned,
             nc::StrCat("link_top_", top_n, "_pinned"), bin_size);

    ++top_n;
    if (top_n == n) {
      break;
    }
  }
}

static std::vector<std::string> PlotLinkUtilizationDeltas() {
  std::chrono::milliseconds bin_size = GetBinSize();

  std::map<StrPair, NumDataVector> data =
      SimpleParseNumericData(FLAGS_input, kLinkUtilizationMetric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);
  std::map<StrPair, NumDataVector> data_pinned =
      SimpleParseNumericData(FLAGS_input_pinned, kLinkUtilizationMetric, ".*",
                             0, std::numeric_limits<uint64_t>::max(), 0);

  std::vector<std::tuple<double, double, std::string>> deltas;
  for (const auto& id_and_rest : data) {
    const StrPair& id = id_and_rest.first;
    const NumDataVector& data_vector = id_and_rest.second;
    if (!nc::ContainsKey(data_pinned, id)) {
      continue;
    }

    const NumDataVector& data_vector_pinned =
        nc::FindOrDieNoPrint(data_pinned, id);

    double avg_utilization =
        GetAverageUtilization(id.second, data_vector, bin_size);
    double avg_utilization_pinned =
        GetAverageUtilization(id.second, data_vector_pinned, bin_size);
    deltas.emplace_back(avg_utilization_pinned, avg_utilization, id.second);
  }

  std::sort(deltas.begin(), deltas.end());

  std::vector<std::pair<double, double>> to_plot;
  std::vector<std::pair<double, double>> to_plot_pinned;

  std::vector<std::string> out;
  for (size_t i = 0; i < deltas.size(); ++i) {
    double avg_utilization;
    double avg_utilization_pinned;
    std::string link;
    std::tie(avg_utilization_pinned, avg_utilization, link) = deltas[i];

    to_plot.emplace_back(avg_utilization, deltas.size() - i);
    to_plot_pinned.emplace_back(avg_utilization_pinned, deltas.size() - i);
    out.emplace_back(link);
  }

  nc::viz::LinePlot scatter_plot(
      {"Link utilization", "utilization", "link (ranked)"});
  scatter_plot.TurnIntoScatterPlot();
  scatter_plot.AddData("regular", to_plot);
  scatter_plot.AddData("pinned", to_plot_pinned);
  scatter_plot.PlotToDir("utilization_deltas");

  return out;
}

static std::vector<std::pair<double, double>> GetSingleDistribution(
    const std::string& input, const std::string& metric) {
  std::map<StrPair, std::vector<nc::DiscreteDistribution<int64_t>>> data =
      SimpleParseDistributionDataNoTimestamps(input, metric, ".*");
  CHECK(data.size() == 1);

  const std::vector<nc::DiscreteDistribution<int64_t>>& data_vector =
      data.begin()->second;
  CHECK(data_vector.size() == 1);

  const nc::DiscreteDistribution<int64_t>& dist = data_vector.front();
  std::vector<int64_t> values = dist.Percentiles(10000);

  std::vector<std::pair<double, double>> out;
  for (size_t i = 0; i < values.size(); ++i) {
    out.emplace_back(values[i], i / 10000.0);
  }

  return out;
}

static void PlotPropDelayDistribution() {
  nc::viz::LinePlot plot(
      {"Propagation delay per aggregate", "delay (ms)", "CDF"});

  plot.AddData("regular", GetSingleDistribution(FLAGS_input, kPropDelayMetric));
  plot.AddData("pinned",
               GetSingleDistribution(FLAGS_input_pinned, kPropDelayMetric));
  plot.PlotToDir("propagation_delays");
}

static void PlotMaxQueueSizeDeltas(const std::vector<std::string>& link_order) {
  std::map<StrPair, NumDataVector> data =
      SimpleParseNumericData(FLAGS_input, kQueueSizeMetric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);
  std::map<StrPair, NumDataVector> data_pinned =
      SimpleParseNumericData(FLAGS_input_pinned, kQueueSizeMetric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);

  std::vector<std::pair<double, double>> to_plot;
  std::vector<std::pair<double, double>> to_plot_pinned;
  for (size_t i = 0; i < link_order.size(); ++i) {
    const std::string& link = link_order[i];

    const NumDataVector& data_vector =
        nc::FindOrDieNoPrint(data, std::make_pair(kQueueSizeMetric, link));
    const NumDataVector& data_vector_pinned = nc::FindOrDieNoPrint(
        data_pinned, std::make_pair(kQueueSizeMetric, link));

    std::chrono::milliseconds queue_size = GetMaxQueueSize(link, data_vector);
    std::chrono::milliseconds queue_size_pinned =
        GetMaxQueueSize(link, data_vector_pinned);

    to_plot.emplace_back(queue_size.count(), link_order.size() - i);
    to_plot_pinned.emplace_back(queue_size_pinned.count(),
                                link_order.size() - i);
  }

  nc::viz::LinePlot scatter_plot(
      {"Queue sizes", "queue size (ms)", "link (ranked by delta utilization)"});
  scatter_plot.TurnIntoScatterPlot();
  scatter_plot.AddData("regular", to_plot);
  scatter_plot.AddData("pinned", to_plot_pinned);
  scatter_plot.PlotToDir("queue_sizes");
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_input.empty()) << "need --input";

  PlotTopNLinks(FLAGS_n);
  PrintUpdateStats();
  PlotPropDelayDistribution();
  if (!FLAGS_input_pinned.empty()) {
    std::vector<std::string> link_order = PlotLinkUtilizationDeltas();
    PlotMaxQueueSizeDeltas(link_order);
  }
  return 0;
}
