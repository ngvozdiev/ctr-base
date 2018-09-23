#ifndef CTR_INFO_H
#define CTR_INFO_H

#include <ncode/net/net_common.h>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <thread>
#include <vector>

#include "../../build/info.pb.h"
#include "../common.h"

namespace ctr {
class RoutingConfiguration;
class TrafficMatrix;
} /* namespace ctr */

namespace ctr {

// Stores topology, traffic matrix and routing information.
class InfoStorage {
 public:
  uint64_t AddTopology(const nc::net::GraphStorage& graph,
                       const std::string& name);

  uint64_t AddTrafficMatrix(uint64_t topology_id, const ctr::TrafficMatrix& tm);

  uint64_t AddRouting(const std::string& routing_system, uint64_t topology_id,
                      uint64_t traffic_matrix_id,
                      const ctr::RoutingConfiguration& routing);

  const info::TopologyInfo* GetTopologyInfoOrNull(uint64_t topology_id) const;

  std::vector<const info::TopologyInfo*> GetAllTopologyInfos() const;

  std::vector<const info::TrafficMatrixInfo*> GetAllTrafficMatrixInfos(
      uint64_t topology_id) const;

  const info::TrafficMatrixInfo* GetTrafficMatrixInfoOrNull(
      uint64_t traffic_matrix_id) const;

  const info::RoutingInfo* GetRoutingInfoOrNull(uint64_t routing_id) const;

  std::vector<const info::RoutingInfo*> GetAllRoutingInfos(
      uint64_t traffic_matrix_id) const;

  const info::RoutingInfo* GetRoutingInfoOrNull(
      uint64_t traffic_matrix_id, const std::string& routing_system) const;

  void DumpToFile(const std::string& file) const;

  void ReadFromFile(const std::string& file);

 private:
  bool IdTaken(uint64_t id) const;

  void PrivateAddTopology(const info::TopologyInfo& topology_info);
  void PrivateAddTrafficMatrix(const info::TrafficMatrixInfo& tm_info);
  void PrivateAddRouting(const info::RoutingInfo& routing_info);

  std::mt19937 rnd_;

  std::map<uint64_t, info::TopologyInfo> topology_infos_;
  std::map<uint64_t, info::TrafficMatrixInfo> traffic_matrix_infos_;
  std::map<uint64_t, info::RoutingInfo> routing_infos_;

  std::map<uint64_t, std::vector<const info::TrafficMatrixInfo*>>
      topology_to_tm_infos_;
  std::map<uint64_t, std::vector<const info::RoutingInfo*>>
      tm_to_routing_infos_;
  std::map<std::string, std::map<uint64_t, const info::RoutingInfo*>>
      routing_system_to_tm_to_routing_infos_;

  std::mutex mu_;
};

// Converts a topology info to a graph.
std::unique_ptr<nc::net::GraphStorage> InfoToGraph(
    const info::TopologyInfo& info);

// Converts a traffic matrix info to a traffic matrix.
std::unique_ptr<ctr::TrafficMatrix> InfoToTM(
    const nc::net::GraphStorage& graph, const info::TrafficMatrixInfo& info);

}  // namespace ctr
#endif
