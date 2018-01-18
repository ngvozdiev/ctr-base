#include <gflags/gflags.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/map_util.h>
#include <ncode/net/net_common.h>
#include <ncode/strutil.h>
#include <ncode/thread_runner.h>
#include <chrono>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "common.h"
#include "demand_matrix_input.h"
#include "opt/ctr.h"
#include "opt/opt.h"
#include "opt/path_provider.h"
#include "topology_input.h"

using namespace std::chrono;

DEFINE_uint64(threads, 4, "Number of parallel threads to run");
DEFINE_bool(run_ctr_link_based, false, "Also run a link-based version of CTR");

struct Input {
  const ctr::DemandMatrixAndFilename* demand_matrix_and_filename;
  const std::vector<std::string>* node_order;
};

namespace ctr {

static constexpr double kLinkScaleStride = 0.02;

static std::vector<double> GetHeadroomVsDelay(const TrafficMatrix& tm) {
  using namespace std::chrono;
  const nc::net::GraphStorage* graph = tm.graph();
  PathProvider path_provider(graph);

  std::vector<double> out;
  for (double link_capacity_multiplier = 1.0; link_capacity_multiplier > 0;
       link_capacity_multiplier -= kLinkScaleStride) {
    if (!tm.ToDemandMatrix()->IsFeasible({}, link_capacity_multiplier)) {
      break;
    }

    CTROptimizer ctr_optimizer(&path_provider, link_capacity_multiplier, false,
                               false);
    std::unique_ptr<RoutingConfiguration> routing = ctr_optimizer.Optimize(tm);
    nc::net::Delay total_delay = routing->TotalPerFlowDelay();
    double value = duration_cast<milliseconds>(total_delay).count() / 1000.0;
    out.emplace_back(value);
  }

  return out;
}

static std::string GetFilename(const std::string& tm_file,
                               const std::string opt_string) {
  std::string tm_base = nc::StringReplace(tm_file, ".demands", "", true);
  return nc::StrCat(tm_base, "_", opt_string, ".rc");
}

static void RecordRoutingConfig(const std::string& out,
                                const std::vector<std::string>& node_order,
                                const RoutingConfiguration& routing) {
  std::string serialized = routing.SerializeToText(node_order);
  LOG(INFO) << "Will write routing config to " << out;
  nc::File::WriteStringToFileOrDie(serialized, out);
}

static void OptAndRecord(const std::string& tm_file,
                         const ctr::TrafficMatrix& tm,
                         const std::vector<std::string>& node_order,
                         const std::string opt_string, Optimizer* opt) {
  std::string out = GetFilename(tm_file, opt_string);
  if (nc::File::Exists(out)) {
    LOG(INFO) << "File exists " << out << ", will skip";
    return;
  }

  auto rc = opt->OptimizeWithTimeAndOptString(tm, opt_string);
  RecordRoutingConfig(out, node_order, *rc);
}

static void RunOptimizers(const Input& input) {
  const std::string& top_file = input.demand_matrix_and_filename->topology_file;
  const std::string& tm_file = input.demand_matrix_and_filename->file;
  nc::lp::DemandMatrix* demand_matrix =
      input.demand_matrix_and_filename->demand_matrix.get();
  const nc::net::GraphStorage* graph = demand_matrix->graph();
  const std::vector<std::string>& node_order = *(input.node_order);
  LOG(ERROR) << "Running " << top_file << " " << tm_file;

  std::unique_ptr<TrafficMatrix> tm =
      TrafficMatrix::DistributeFromDemandMatrix(*demand_matrix);
  if (!nc::ContainsKey(demand_matrix->properties(), "headroom_vs_delay")) {
    std::vector<double> headroom_vs_delay = GetHeadroomVsDelay(*tm);
    std::string to_record = nc::Join(headroom_vs_delay, ",");

    demand_matrix->UpdateProperty("headroom_vs_delay", to_record);
    demand_matrix->UpdateProperty("headroom_vs_delay_step",
                                  nc::StrCat(kLinkScaleStride));
    LOG(INFO) << "Will overwrite " << tm_file << " with headroom vs delay info";
    demand_matrix->ToRepetitaFileOrDie(node_order, tm_file);
  }

  PathProvider path_provider(graph);
  CTROptimizer ctr_optimizer(&path_provider, 1.0, false, false);
  CTROptimizer ctr_optimizer_no_flow_counts(&path_provider, 1.0, false, true);
  MinMaxOptimizer minmax_low_delay_optimizer(&path_provider, 1.0, true);
  MinMaxPathBasedOptimizer minmax_ksp_optimizer(&path_provider, 1.0, true, 10);
  B4Optimizer b4_optimizer(&path_provider, false, 1.0);
  CTRLinkBased ctr_link_based_optimizer(&path_provider, 1.0);
  ShortestPathOptimizer sp_optimizer(&path_provider);

  OptAndRecord(tm_file, *tm, node_order, "CTRNOCACHE", &ctr_optimizer);
  OptAndRecord(tm_file, *tm, node_order, "CTR", &ctr_optimizer);
  OptAndRecord(tm_file, *tm, node_order, "CTRNFC",
               &ctr_optimizer_no_flow_counts);
  OptAndRecord(tm_file, *tm, node_order, "MinMaxLD",
               &minmax_low_delay_optimizer);
  OptAndRecord(tm_file, *tm, node_order, "MinMaxK10", &minmax_ksp_optimizer);
  OptAndRecord(tm_file, *tm, node_order, "B4", &b4_optimizer);
  OptAndRecord(tm_file, *tm, node_order, "SP", &sp_optimizer);
  if (FLAGS_run_ctr_link_based) {
    OptAndRecord(tm_file, *tm, node_order, "CTRLB", &ctr_link_based_optimizer);
  }
}

}  // namespace ctr

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::tie(topologies, matrices) = ctr::GetDemandMatrixInputs(true);

  std::map<std::string, const ctr::TopologyAndFilename*> topologies_by_name;
  for (const auto& topology : topologies) {
    topologies_by_name[topology.file] = &topology;
  }

  std::vector<Input> to_process;
  for (const auto& matrix : matrices) {
    const std::string& top_file = matrix.topology_file;
    const ctr::TopologyAndFilename* top_ptr =
        nc::FindOrDie(topologies_by_name, top_file);
    to_process.push_back({&matrix, &top_ptr->node_order});
  }

  nc::RunInParallel<Input>(to_process, [](const Input& input) {
    ctr::RunOptimizers(input);
  }, FLAGS_threads);
}
