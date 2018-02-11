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
DEFINE_uint64(line_plot_bin_size, 1,
              "Will bin every N values based on "
              "timestamp and plot the means in each bin.");
DEFINE_uint64(n, 5, "Will plot the top N links");

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

static nc::net::Bandwidth GetLinkRate(const StrPair& link_src_and_dst) {
  std::map<StrPair, std::vector<double>> data =
      SimpleParseNumericDataNoTimestamps(FLAGS_input, kLinkRateMetric, ".*");
  const std::vector<double>& values =
      nc::FindOrDieNoPrint(data, link_src_and_dst);
  CHECK(values.size() == 1);
  return nc::net::Bandwidth::FromMBitsPerSecond(values[0]);
}

static double GetAverageUtilization(const StrPair& link_src_and_dst,
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

static void PlotLink(const StrPair& link_src_and_dst,
                     const NumDataVector& link_utilization,
                     const std::string& output,
                     std::chrono::milliseconds bin_size) {
  nc::net::Bandwidth link_rate = GetLinkRate(link_src_and_dst);

  viz::PlotParameters2D plot_params;
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
  for (auto it = data.cbegin(); it != data.cend(); ++it) {
    const StrPair& src_and_dst = it->first;
    const NumDataVector& data_vector = it->second;

    double avg_utilization =
        GetAverageUtilization(src_and_dst, data_vector, bin_size);
    top_links.emplace(avg_utilization, it);
  }

  uint32_t top_n = 0;
  for (auto rit = top_links.rbegin(); rit != top_links.rend(); ++rit) {
    double avg_utilization = rit->first;
    auto data_it = rit->second;
    const StrPair& src_and_dst = data_it->first;
    const NumDataVector& data_vector = data_it->second;

    std::string out = nc::StrCat("link_top_", top_n);
    PlotLink(src_and_dst, data_vector, out, bin_size);

    ++top_n;
    if (top_n == n) {
      break;
    }
  }
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_input.empty()) << "need --input";

  PlotTopNLinks(FLAGS_n);
  return 0;
}
