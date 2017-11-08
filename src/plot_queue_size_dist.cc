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

// Returns for each link (queue) in the topology its service rate.
static std::map<std::string, double> GetRatesMbps() {
  std::map<std::pair<std::string, std::string>, NumDataVector> data =
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
  std::map<std::pair<std::string, std::string>, NumDataVector> data =
      SimpleParseNumericDataNoTimestamps(FLAGS_input, kRecordPeriodMetric,
                                         ".*");
  CHECK(data.size() == 1);
  const NumDataVector& data_vector = data.begin()->second;
  CHECK(data_vector.size() == 1);
  return std::chrono::milliseconds(data_vector.front());
}

static std::vector<double> LinkUtilization(const std::string& id,
                                           double rate_Mbps) {}

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
    LOG(INFO) << "P " << i / static_cast<double>(values.size()) << " "
              << values[i];
    data_series.data.emplace_back(values[i], i);
  }

  nc::viz::PythonGrapher grapher("queue_size_plot_out");
  grapher.PlotLine({}, {data_series});
  return 0;
}
