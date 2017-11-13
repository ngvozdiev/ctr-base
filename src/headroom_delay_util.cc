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

// Produces a plot of increase in headroom vs increase in total propagation
// delay for a topology/tm combo.
static void HeadroomVsDelayEval(const ctr::OptEvalInput& input, bool verbose) {
  using namespace ctr;
  using namespace std::chrono;
  LOG(INFO) << "Processing topology " << input.topology_file << " tm "
            << input.tm_file;

  const nc::lp::DemandMatrix& demand_matrix = *(input.demand_matrix);
  const nc::net::GraphStorage* graph = demand_matrix.graph();

  const std::string& top_file = input.topology_file;
  const std::string& tm_file = input.tm_file;

  auto tm = TrafficMatrix::ProportionalFromDemandMatrix(demand_matrix);
  std::vector<double> values;
  for (double link_multiplier = 1.0; link_multiplier > 0.0;
       link_multiplier -= FLAGS_stride_size) {
    PathProvider path_provider(graph);
    CTROptimizer ctr_optimizer(&path_provider, link_multiplier, false, false);
    std::unique_ptr<RoutingConfiguration> routing = ctr_optimizer.Optimize(*tm);

    double max_utilization = routing->MaxLinkUtilization();
    LOG_IF(INFO, verbose) << "link multiplier " << link_multiplier
                          << " max utilization " << max_utilization
                          << " scale factor "
                          << demand_matrix.MaxCommodityScaleFractor(
                                 link_multiplier);

    if (max_utilization > link_multiplier + 0.001) {
      break;
    }

    double value = duration_cast<seconds>(routing->TotalPerFlowDelay()).count();
    values.emplace_back(value);
  }

  // The lowest delay is achieved at link_multiplier=1.0. Will normalize
  // everything by that.
  CHECK(!values.empty());
  double lowest_delay = values.front();
  auto* metric_handle = total_delay_fraction->GetHandle(top_file, tm_file);
  for (auto value : values) {
    metric_handle->AddValue(value / lowest_delay);
  }
}

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
