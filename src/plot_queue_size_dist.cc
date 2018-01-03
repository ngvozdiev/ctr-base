#include <gflags/gflags.h>
#include <stddef.h>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ncode/logging.h"
#include "ncode/stats.h"
#include "ncode/viz/grapher.h"
#include "metrics/metrics_parser.h"

DEFINE_string(input, "", "The metrics file.");
DEFINE_uint64(bin_size_sec, 10,
              "Over what period of time to average link utilization.");

static constexpr char kQueueSizeMetric[] = "queue_size_ms";
static constexpr char kPropDelayMetric[] = "propagation_delay_ms";

using namespace nc;
using namespace nc::metrics::parser;

using UintDataVector = std::vector<nc::DiscreteDistribution<int64_t>>;
using NumDataVector = std::vector<double>;
using StrPair = std::pair<std::string, std::string>;

static void PlotDistribution(const std::string& metric,
                             const std::string& out) {
  std::map<StrPair, UintDataVector> data =
      SimpleParseDistributionDataNoTimestamps(FLAGS_input, metric, ".*");
  CHECK(data.size() == 1);

  const UintDataVector& data_vector = data.begin()->second;
  CHECK(data_vector.size() == 1);

  const nc::DiscreteDistribution<int64_t>& dist = data_vector.front();
  std::vector<int64_t> values = dist.Percentiles(10000);

  nc::viz::DataSeries2D data_series;
  for (size_t i = 0; i < values.size(); ++i) {
    data_series.data.emplace_back(values[i], i);
  }

  nc::viz::LinePlot plot;
  plot.AddData(data_series);
  plot.PlotToDir(out);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_input.empty()) << "need --input";

  PlotDistribution(kQueueSizeMetric, "queue_size_plot_out");
  PlotDistribution(kPropDelayMetric, "propagation_delay_plot_out");
  return 0;
}
