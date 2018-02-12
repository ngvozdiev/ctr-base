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
DEFINE_uint64(n, 10, "Will plot the top N links");

static constexpr char kLinkUtilizationMetric[] = "link_utilization";
static constexpr char kLinkRateMetric[] = "link_rate_Mbps";
static constexpr char kBinSizeMetric[] = "bin_size_ms";

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

  std::map<StrPair, NumDataVector> data =
      SimpleParseNumericData(FLAGS_input, kLinkUtilizationMetric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);
  std::map<double, std::map<StrPair, NumDataVector>::const_iterator> top_links;
  std::vector<double> all_utilizations;
  for (auto it = data.cbegin(); it != data.cend(); ++it) {
    const std::string& src_and_dst = it->first.second;
    const NumDataVector& data_vector = it->second;

    double avg_utilization =
        GetAverageUtilization(src_and_dst, data_vector, bin_size);
    all_utilizations.emplace_back(avg_utilization);
    top_links.emplace(avg_utilization, it);
  }

  uint32_t top_n = 0;
  for (auto rit = top_links.rbegin(); rit != top_links.rend(); ++rit) {
    double avg_utilization = rit->first;
    auto data_it = rit->second;
    const std::string& src_and_dst = data_it->first.second;
    const NumDataVector& data_vector = data_it->second;

    std::string out = nc::StrCat("link_top_", top_n);

    LOG(INFO) << src_and_dst << " " << avg_utilization;
    PlotLink(src_and_dst, data_vector, out, bin_size);

    ++top_n;
    if (top_n == n) {
      break;
    }
  }

  nc::viz::CDFPlot utilization_cdf;
  utilization_cdf.AddData("", all_utilizations);
  utilization_cdf.PlotToDir("all_utilizations");
}

static void PlotLinkUtilizationDeltas() {
  std::chrono::milliseconds bin_size = GetBinSize();

  std::map<StrPair, NumDataVector> data =
      SimpleParseNumericData(FLAGS_input, kLinkUtilizationMetric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);
  std::map<StrPair, NumDataVector> data_pinned =
      SimpleParseNumericData(FLAGS_input_pinned, kLinkUtilizationMetric, ".*",
                             0, std::numeric_limits<uint64_t>::max(), 0);

  std::vector<std::pair<double, bool>> deltas;
  for (const auto& id_and_rest : data) {
    const StrPair& id = id_and_rest.first;
    const NumDataVector& data_vector = id_and_rest.second;
    const NumDataVector& data_vector_pinned =
        nc::FindOrDieNoPrint(data_pinned, id);

    double avg_utilization =
        GetAverageUtilization(id.second, data_vector, bin_size);
    double avg_utilization_pinned =
        GetAverageUtilization(id.second, data_vector_pinned, bin_size);
    double delta = avg_utilization - avg_utilization_pinned;
    bool top_utilization = avg_utilization_pinned > 0.999;
    deltas.emplace_back(delta, top_utilization);
  }

  std::sort(deltas.begin(), deltas.end());

  std::vector<std::pair<double, double>> to_plot;
  std::vector<std::pair<double, double>> to_plot_top_utilization;
  for (size_t i = 0; i < deltas.size(); ++i) {
    double delta;
    bool top_utilization;
    std::tie(delta, top_utilization) = deltas[i];

    if (top_utilization) {
      to_plot_top_utilization.emplace_back(delta, i);
    } else {
      to_plot.emplace_back(delta, i);
    }
  }

  nc::viz::LinePlot scatter_plot(
      {"Delta link utilization", "utilization - utilization pinned", "link"});
  scatter_plot.TurnIntoScatterPlot();
  scatter_plot.AddData("other", to_plot);
  scatter_plot.AddData("top_utilization", to_plot_top_utilization);
  scatter_plot.PlotToDir("utilization_deltas");
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_input.empty()) << "need --input";

  PlotTopNLinks(FLAGS_n);
  if (!FLAGS_input_pinned.empty()) {
    PlotLinkUtilizationDeltas();
  }
  return 0;
}
