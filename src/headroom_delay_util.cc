#include <gflags/gflags.h>
#include <chrono>
#include <iostream>
#include <memory>
#include <tuple>
#include <type_traits>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/thread_runner.h"
#include "common.h"
#include "metrics/metrics.h"
#include "opt/ctr.h"
#include "opt/path_provider.h"
#include "opt_eval.h"

using namespace std::chrono;

DEFINE_uint64(threads, 4, "Number of parallel threads to run");
DEFINE_double(stride_size, 0.2, "How much to increase headroom each iteration");

static auto* stride_size_metric =
    nc::metrics::DefaultMetricManager() -> GetThreadSafeMetric<double>(
        "stride_size", "How much further from 1.0 each step is");

static auto* total_delay_fraction =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string>(
            "total_delay_fraction", "Fraction of total delay at no headroom",
            "Topology", "Traffic matrix");



using TopologyAndMatrix = std::tuple<std::string, std::string, double>;
int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::metrics::InitMetrics();
  stride_size_metric->GetHandle()->AddValue(FLAGS_stride_size);

  // Will switch off timestamps.
  auto timestamp_provider =
      ::nc::make_unique<nc::metrics::NullTimestampProvider>();
  nc::metrics::DefaultMetricManager()->set_timestamp_provider(
      std::move(timestamp_provider));

  std::vector<std::unique_ptr<nc::net::GraphStorage>> graphs;
  std::vector<ctr::OptEvalInput> to_process;
  std::tie(graphs, to_process) = ctr::GetOptEvalInputs();

  nc::RunInParallel<ctr::OptEvalInput>(
      to_process, [&to_process](const ctr::OptEvalInput& input) {
        HeadroomVsDelayEval(input, to_process.size() == 1);
      }, FLAGS_threads);
}
