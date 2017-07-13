#include <gflags/gflags.h>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "common.h"
#include "opt/ctr.h"
#include "opt/opt.h"
#include "opt/opt_compare.h"
#include "opt/path_provider.h"

DEFINE_string(topology, "", "A file");
DEFINE_string(traffic_matrix, "", "A traffic matrix file");
DEFINE_double(link_capacity_scale, 1.0, "By how much to scale all links");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_string(output, "opt_eval_out", "Output directory");

namespace ctr {

std::unique_ptr<TrafficMatrix> FromDemandMatrix(
    const nc::lp::DemandMatrix& demand_matrix) {
  std::map<AggregateId, DemandAndFlowCount> demands_and_counts;
  for (const auto& element : demand_matrix.elements()) {
    AggregateId id(element.src, element.dst);
    demands_and_counts[id] = {element.demand, 1000ul};
  }

  return nc::make_unique<TrafficMatrix>(demand_matrix.graph(),
                                        demands_and_counts);
}

void RunOptimizer(Optimizer* optimizer, RoutingConfigInfo* out,
                  const nc::lp::DemandMatrix& demand_matrix) {
  std::unique_ptr<RoutingConfiguration> routing =
      optimizer->Optimize(*FromDemandMatrix(demand_matrix));
  std::cout << routing->ToString();
  out->Add(*routing);
}

}  // namespace ctr

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_traffic_matrix.empty());
  CHECK(!FLAGS_topology.empty());

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology), &node_order);
  builder.RemoveMultipleLinks();
  builder.ScaleCapacity(FLAGS_link_capacity_scale);
  nc::net::GraphStorage graph(builder);

  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
      nc::lp::DemandMatrix::LoadRepetitaOrDie(
          nc::File::ReadFileToStringOrDie(FLAGS_traffic_matrix), node_order,
          &graph);
  demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);

  ctr::PathProvider path_provider(&graph);
  ctr::CTROptimizer ctr_optimizer(&path_provider, 1.0, true);
  ctr::B4Optimizer b4_optimizer(&path_provider, false);
  ctr::MinMaxOptimizer minmax_optimizer(&path_provider);

  ctr::RoutingConfigInfo ctr_info;
  ctr::RoutingConfigInfo b4_info;
  ctr::RoutingConfigInfo minmax_info;

  ctr::RunOptimizer(&ctr_optimizer, &ctr_info, *demand_matrix);
  ctr::RunOptimizer(&b4_optimizer, &b4_info, *demand_matrix);
  ctr::RunOptimizer(&minmax_optimizer, &minmax_info, *demand_matrix);

  ctr_info.AddPrefixToFieldNames("ctr_");
  b4_info.AddPrefixToFieldNames("b4_");
  minmax_info.AddPrefixToFieldNames("minmax_");

  ctr_info.Combine(b4_info).Combine(minmax_info);
  ctr_info.ToDisk(FLAGS_output, true);
}
