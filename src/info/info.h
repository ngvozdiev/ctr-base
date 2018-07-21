#ifndef CTR_INFO_H
#define CTR_INFO_H

#include <ncode/net/net_common.h>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include "../../build/info.pb.h"

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

  uint64_t AddRouting(uint64_t topology_id, uint64_t traffic_matrix_id,
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
  uint64_t GetId();

  std::mt19937 rnd_;
  std::uniform_int_distribution<uint64_t> distribution_;

  std::map<uint64_t, info::TopologyInfo> topology_infos_;
  std::map<uint64_t, info::TrafficMatrixInfo> traffic_matrix_infos_;
  std::map<uint64_t, info::RoutingInfo> routing_infos_;
};

}  // namespace ctr
#endif
