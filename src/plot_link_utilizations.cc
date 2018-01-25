#include <gflags/gflags.h>
#include <ncode/common.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/net/net_common.h>
#include <ncode/strutil.h>
#include <ncode/viz/grapher.h>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "common.h"
#include "demand_matrix_input.h"
#include "opt/path_provider.h"
#include "topology_input.h"

DEFINE_string(opt, "CTR,MinMaxLD", "The optimizers to plot");

static std::string GetFilename(const std::string& tm_file,
                               const std::string opt_string) {
  std::string tm_base = nc::StringReplace(tm_file, ".demands", "", true);
  return nc::StrCat(tm_base, "_", opt_string, ".rc");
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
  const std::vector<std::string>& node_order = topologies[0].node_order;
  const std::string& demand_file = matrices[0].file;
  const std::string& topology_file = topologies[0].file;

  nc::viz::CDFPlot plot({nc::StrCat("Link utilizations for ",
                                    nc::File::ExtractFileName(demand_file)),
                         "link utilization"});
  std::vector<std::string> opts = nc::Split(FLAGS_opt, ",", true);
  for (const std::string& opt : opts) {
    ctr::PathProvider path_provider(graph);
    std::string rc_filename = GetFilename(demand_file, opt);
    if (!nc::File::Exists(rc_filename)) {
      LOG(INFO) << "Missing " << rc_filename;
      continue;
    }

    std::string rc_serialized = nc::File::ReadFileToStringOrDie(rc_filename);
    std::unique_ptr<ctr::TrafficMatrix> tm =
        ctr::TrafficMatrix::DistributeFromDemandMatrix(*demand_matrix);
    LOG(INFO) << "Will parse " << demand_file << " at " << topology_file
              << " opt " << opt;
    auto rc = ctr::RoutingConfiguration::LoadFromSerializedText(
        *tm, node_order, rc_serialized, &path_provider);

    nc::net::GraphLinkMap<double> utilizations = rc->LinkUtilizations();
    std::vector<double> to_plot;
    for (const auto& link_and_util : utilizations) {
      to_plot.emplace_back(*link_and_util.second);
    }

    double total = std::accumulate(to_plot.begin(), to_plot.end(), 0.0);
    double mean = total / to_plot.size();
    std::string label =
        nc::StrCat(opt, " (mean ", nc::ToStringMaxDecimals(mean, 2), ")");
    plot.AddData(label, to_plot);
  }

  plot.PlotToDir("utilization");
}
