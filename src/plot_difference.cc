#include <gflags/gflags.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/net/net_common.h>
#include <ncode/strutil.h>
#include <ncode/viz/web_page.h>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "common.h"
#include "demand_matrix_input.h"
#include "opt/path_provider.h"
#include "topology_input.h"

DEFINE_string(from, "CTR", "Difference from");
DEFINE_string(to, "B4", "Difference to");

static std::string GetFilename(const std::string& tm_file,
                               const std::string opt_string) {
  std::string tm_base = nc::StringReplace(tm_file, ".demands", "", true);
  return nc::StrCat(tm_base, "_", opt_string, ".rc");
}

static std::unique_ptr<ctr::RoutingConfiguration> GetRC(
    const ctr::TopologyAndFilename& top,
    const ctr::DemandMatrixAndFilename& demand, const std::string& opt,
    ctr::PathProvider* path_provider) {
  const std::string& demand_file = demand.file;
  const nc::lp::DemandMatrix& demand_matrix = *(demand.demand_matrix);

  std::string rc_filename = GetFilename(demand_file, opt);
  if (!nc::File::Exists(rc_filename)) {
    LOG(FATAL) << "Missing " << rc_filename;
  }

  std::string rc_serialized = nc::File::ReadFileToStringOrDie(rc_filename);
  std::unique_ptr<ctr::TrafficMatrix> tm =
      ctr::TrafficMatrix::DistributeFromDemandMatrix(demand_matrix);
  LOG(INFO) << "Will parse " << demand_file << " opt " << opt;

  const std::vector<std::string>& node_order = top.node_order;
  return ctr::RoutingConfiguration::LoadFromSerializedText(
      *tm, node_order, rc_serialized, path_provider);
}

static void OverloadedLinks(const ctr::RoutingConfiguration& rc) {
  const nc::net::GraphStorage* graph = rc.graph();
  for (const auto& link_and_util : rc.LinkUtilizations()) {
    if (*link_and_util.second > 1) {
      LOG(INFO) << "Link "
                << graph->GetLink(link_and_util.first)->ToStringNoPorts()
                << " load " << *(link_and_util.second);
    }
  }
}

static void DumpMaxKAggregate(const ctr::RoutingConfiguration& rc) {
  nc::net::Delay max_abs_delay = nc::net::Delay::zero();
  const ctr::AggregateId* max_abs_delay_id = nullptr;
  double max_rel_delay = 0.0;
  const ctr::AggregateId* max_rel_delay_id = nullptr;
  double max_hc = std::numeric_limits<double>::min();
  const ctr::AggregateId* max_hc_id = nullptr;

  const nc::net::GraphStorage* graph = rc.graph();
  const std::map<ctr::AggregateId, std::vector<ctr::RouteAndFraction>>&
      aggregates = rc.routes();
  for (const auto& id_and_routes : aggregates) {
    const ctr::AggregateId& id = id_and_routes.first;
    const std::vector<ctr::RouteAndFraction>& routes = id_and_routes.second;

    nc::net::Delay sp_delay = id.GetSPDelay(*graph);
    double sp_hc = id.GetSP(*graph)->links().size();
    for (const auto& route : routes) {
      nc::net::Delay delay = route.first->delay();
      nc::net::Delay abs_delay = delay - sp_delay;
      double hc = route.first->links().size();

      if (abs_delay > max_abs_delay) {
        max_abs_delay = abs_delay;
        max_abs_delay_id = &id;
      }

      double rel_delay = static_cast<double>(delay.count()) / sp_delay.count();
      if (rel_delay > max_rel_delay) {
        max_rel_delay = rel_delay;
        max_rel_delay_id = &id;
      }

      double hc_delta = hc - sp_hc;
      if (sp_hc == 1.0 && hc_delta > max_hc) {
        max_hc = hc_delta;
        max_hc_id = &id;
      }
    }
  }

  CHECK(max_abs_delay_id != nullptr);
  CHECK(max_rel_delay_id != nullptr);
  CHECK(max_hc_id != nullptr);

  LOG(INFO) << "Max abs "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   max_abs_delay)
                   .count()
            << " aggregate " << max_abs_delay_id->ToString(*rc.graph());
  LOG(INFO) << "Max rel " << max_rel_delay << " aggregate "
            << max_rel_delay_id->ToString(*rc.graph());
  LOG(INFO) << "Max hc " << max_hc << " aggregate "
            << max_hc_id->ToString(*rc.graph());
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::tie(topologies, matrices) = ctr::GetDemandMatrixInputs(false);
  CHECK(topologies.size() == 1) << "Only one topology should match";
  CHECK(matrices.size() == 1) << "Only one TM should match";

  const nc::lp::DemandMatrix* demand_matrix = matrices[0].demand_matrix.get();
  const nc::net::GraphStorage* graph = demand_matrix->graph();

  ctr::PathProvider path_provider(graph);
  auto from_rc = GetRC(topologies[0], matrices[0], FLAGS_from, &path_provider);
  auto to_rc = GetRC(topologies[0], matrices[0], FLAGS_to, &path_provider);
  DumpMaxKAggregate(*to_rc);

  ctr::RoutingConfigurationDelta delta = from_rc->GetDifference(*to_rc);
  LOG(INFO) << delta.ToString(*graph);
  LOG(INFO) << "overloaded from " << from_rc->OverloadedAggregates().size()
            << " to " << to_rc->OverloadedAggregates().size();
  //
  //  nc::viz::HtmlPage from_out;
  //  from_rc->ToHTML(&from_out);
  //  nc::File::WriteStringToFileOrDie(from_out.Construct(), "from_out.html");

  nc::viz::HtmlPage to_out;
  to_rc->ToHTML(&to_out);
  nc::File::WriteStringToFileOrDie(to_out.Construct(), "to_out.html");
}
