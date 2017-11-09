#include <gflags/gflags.h>
#include <stddef.h>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ncode_common/src/logging.h"
#include "ncode_common/src/stats.h"
#include "ncode_common/src/viz/grapher.h"
#include "metrics/metrics_parser.h"

DEFINE_string(input, "", "The metrics file.");
DEFINE_uint64(bin_size_sec, 10,
              "Over what period of time to average link utilization.");

static constexpr char kQueueSizeMetric[] = "queue_size_ms";
static constexpr char kBytesSeenMetric[] = "net_queue_bytes_seen";
static constexpr char kQueueRateMetric[] = "queue_rate_Mbps";
static constexpr char kRecordPeriodMetric[] = "record_period";

using namespace nc;
using namespace nc::metrics::parser;

using UintDataVector = std::vector<nc::DiscreteDistribution<int64_t>>;
using NumDataVector = std::vector<double>;
using StrPair = std::pair<std::string, std::string>;

// Returns for each link (queue) in the topology its service rate.
static std::map<std::string, double> GetRatesMbps() {
  std::map<StrPair, NumDataVector> data =
      SimpleParseNumericDataNoTimestamps(FLAGS_input, kQueueRateMetric, ".*");

  std::map<std::string, double> out;
  for (const auto& id_and_data : data) {
    const std::string& metric_id = id_and_data.first.second;
    CHECK(id_and_data.second.size() == 1);
    out[metric_id] = id_and_data.second.front();
  }

  return out;
}

static std::chrono::milliseconds GetRecordPeriod() {
  std::map<StrPair, NumDataVector> data = SimpleParseNumericDataNoTimestamps(
      FLAGS_input, kRecordPeriodMetric, ".*");
  CHECK(data.size() == 1);
  const NumDataVector& data_vector = data.begin()->second;
  CHECK(data_vector.size() == 1);
  return std::chrono::milliseconds(static_cast<uint64_t>(data_vector.front()));
}

// Returns n-1 elements with the differences between values.
static std::vector<double> Diff(const std::vector<double>& v) {
  std::vector<double> out(v.size() - 1);
  for (size_t i = 0; i < v.size() - 1; ++i) {
    double delta = v[i + 1] - v[i];
    out[i] = delta;
  }

  return out;
}

// Returns the link utilization in time.
static std::vector<double> LinkUtilization(
    const std::vector<double>& bytes_seen,
    std::chrono::milliseconds record_period, double rate_Mbps) {
  double num_periods_per_second = 1000.0 / record_period.count();
  size_t num_records_per_bin = FLAGS_bin_size_sec * num_periods_per_second;

  // Each bin will contain the mean value of all readings within the bin.
  std::vector<double> binned(bytes_seen.begin(), bytes_seen.end());
  nc::Bin(num_records_per_bin, &binned);

  for (auto& v : binned) {
    v *= 8;                       // Bytes to bits.
    v *= num_periods_per_second;  // Per period to per second.
    v /= 1000000.0;               // bits per second to Mbps.
    v /= rate_Mbps;  // rate to utilization.
    CHECK(v >= 0);
    CHECK(v <= 1.0);
  }

  //  LOG(FATAL) << "Bin count " << binned.size() << " num records per bin "
  //             << num_records_per_bin << " record period "
  //             << record_period.count() << "ms";

  return binned;
}

// For each link returns a vector with 101 values---the 100 percentiles in the
// distribution of utilization.
static std::map<std::string, std::vector<double>> AllLinkUtilizations() {
  std::chrono::milliseconds record_period = GetRecordPeriod();
  std::map<StrPair, NumDataVector> bytes_seen_data =
      SimpleParseNumericDataNoTimestamps(FLAGS_input, kBytesSeenMetric, ".*");
  std::map<std::string, double> rates = GetRatesMbps();

  std::map<std::string, std::vector<double>> out;
  for (const auto& id_and_data : bytes_seen_data) {
    const std::string& link_id = id_and_data.first.second;
    const std::vector<double>& bytes_seen = id_and_data.second;
    double rate_Mbps = nc::FindOrDie(rates, link_id);

    std::vector<double> link_utilization =
        LinkUtilization(Diff(bytes_seen), record_period, rate_Mbps);
    out.emplace(link_id, nc::Percentiles(&link_utilization, 100));
  }

  return out;
}

static std::vector<double> GetLinkUtilizationsAtPercentile(
    const std::map<std::string, std::vector<double>>& utilizations, size_t p) {
  std::vector<double> out;
  for (const auto& link_and_utilization : utilizations) {
    std::vector<double> utilization_p = link_and_utilization.second;
    out.emplace_back(utilization_p[p]);
  }

  return out;
}

static void PlotLinkUtilizations() {
  std::map<std::string, std::vector<double>> link_utilizations =
      AllLinkUtilizations();

  nc::viz::PythonGrapher grapher("link_utilization_plot_out");
  grapher.PlotCDF(
      {}, {{"max", GetLinkUtilizationsAtPercentile(link_utilizations, 100)},
           {"median", GetLinkUtilizationsAtPercentile(link_utilizations, 50)}});
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_input.empty()) << "need --input";

  std::map<std::pair<std::string, std::string>, UintDataVector> sp_data =
      SimpleParseDistributionDataNoTimestamps(FLAGS_input, kQueueSizeMetric,
                                              ".*");
  CHECK(sp_data.size() == 1);

  const UintDataVector& data_vector = sp_data.begin()->second;
  CHECK(data_vector.size() == 1);

  const nc::DiscreteDistribution<int64_t>& dist = data_vector.front();
  std::vector<int64_t> values = dist.Percentiles(10000);

  nc::viz::DataSeries2D data_series;
  for (size_t i = 0; i < values.size(); ++i) {
    data_series.data.emplace_back(values[i], i);
  }

  nc::viz::PythonGrapher grapher("queue_size_plot_out");
  grapher.PlotLine({}, {data_series});

  PlotLinkUtilizations();
  return 0;
}
