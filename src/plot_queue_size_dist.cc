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

static constexpr char kQueueSizeMetric[] = "queue_size_ms";

using namespace nc;
using namespace nc::metrics::parser;

using DataVector = std::vector<nc::DiscreteDistribution<int64_t>>;

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_input.empty()) << "need --input";

  std::map<std::pair<std::string, std::string>, DataVector> sp_data =
      SimpleParseDistributionDataNoTimestamps(FLAGS_input, kQueueSizeMetric,
                                              ".*");
  CHECK(sp_data.size() == 1);

  const DataVector& data_vector = sp_data.begin()->second;
  CHECK(data_vector.size() == 1);

  const nc::DiscreteDistribution<int64_t>& dist = data_vector.front();
  std::vector<int64_t> values = dist.Percentiles(10000);

  nc::viz::DataSeries1D data_series;
  for (size_t i = 0; i < values.size(); ++i) {
    data_series.data.emplace_back(values[i]);
  }

  nc::viz::PythonGrapher grapher("queue_size_plot_out");
  grapher.PlotCDF({}, {data_series});
  return 0;
}
