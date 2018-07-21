#include <gflags/gflags.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/map_util.h>
#include <ncode/net/net_common.h>
#include <ncode/strutil.h>
#include <stdint.h>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "../common.h"
#include "../demand_matrix_input.h"
#include "../opt/path_provider.h"
#include "../topology_input.h"
#include "info.h"

DEFINE_string(opt, "SP,B4,MinMaxK10,CTR,MinMaxLD", "The optimizers");
DEFINE_string(output, "info.pb", "All infos");

static std::string GetFilename(const std::string& tm_file,
                               const std::string opt_string) {
  std::string tm_base = nc::StringReplace(tm_file, ".demands", "", true);
  return nc::StrCat(tm_base, "_", opt_string, ".rc");
}

static std::string GetTopologyNameFromFilename(const std::string& file) {
  std::string filename = nc::File::ExtractFileName(file);
  std::string::size_type last_dot = filename.find_last_of('.');
  return filename.substr(0, last_dot);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::vector<std::string> opts = nc::Split(FLAGS_opt, ",", true);

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::tie(topologies, matrices) = ctr::GetDemandMatrixInputs(false);

  std::map<const nc::net::GraphStorage*, std::vector<std::string>>
      topology_to_node_order;
  std::map<const nc::net::GraphStorage*, std::string> topology_to_name;
  for (const auto& topology_and_filename : topologies) {
    const nc::net::GraphStorage* graph = topology_and_filename.graph.get();
    topology_to_name[graph] =
        GetTopologyNameFromFilename(topology_and_filename.file);
    topology_to_node_order[graph] = topology_and_filename.node_order;
  }

  std::map<const nc::net::GraphStorage*,
           std::vector<const nc::lp::DemandMatrix*>> data_grouped;
  std::map<const nc::lp::DemandMatrix*, std::string> demand_to_filename;
  for (const auto& demand_and_filename : matrices) {
    const nc::lp::DemandMatrix* demand_matrix =
        demand_and_filename.demand_matrix.get();
    demand_to_filename[demand_matrix] = demand_and_filename.file;
    data_grouped[demand_matrix->graph()].emplace_back(demand_matrix);
  }

  CHECK(!FLAGS_output.empty()) << "Need output file";

  ctr::InfoStorage info_storage;
  for (const auto& graph_and_demands : data_grouped) {
    const nc::net::GraphStorage* graph = graph_and_demands.first;
    const std::vector<std::string>& node_order =
        nc::FindOrDie(topology_to_node_order, graph);
    ctr::PathProvider path_provider(graph);

    uint64_t topology_id = info_storage.AddTopology(
        *graph, nc::FindOrDie(topology_to_name, graph));
    for (const nc::lp::DemandMatrix* demand : graph_and_demands.second) {
      const std::string& filename = nc::FindOrDie(demand_to_filename, demand);
      LOG(INFO) << "Processing " << filename;

      std::unique_ptr<ctr::TrafficMatrix> tm =
          ctr::TrafficMatrix::DistributeFromDemandMatrix(*demand);
      uint64_t traffic_matrix_id =
          info_storage.AddTrafficMatrix(topology_id, *tm);

      for (const std::string& opt : opts) {
        std::string rc_filename = GetFilename(filename, opt);
        if (!nc::File::Exists(rc_filename)) {
          LOG(INFO) << "Missing " << rc_filename;
          continue;
        }

        std::string rc_serialized =
            nc::File::ReadFileToStringOrDie(rc_filename);
        auto rc = ctr::RoutingConfiguration::LoadFromSerializedText(
            *tm, node_order, rc_serialized, &path_provider);
        info_storage.AddRouting(topology_id, traffic_matrix_id, *rc);
      }
    }
  }
  info_storage.DumpToFile(FLAGS_output);
}
