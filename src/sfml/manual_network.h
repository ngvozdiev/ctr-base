#ifndef MANUAL_NET_H
#define MANUAL_NET_H

#include <chrono>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/htsim/network.h"
#include "ncode_common/src/htsim/packet.h"
#include "ncode_common/src/htsim/queue.h"
#include "ncode_common/src/net/net_common.h"
#include "../common.h"

namespace ctr {
namespace controller {
class DeviceFactory;
} /* namespace controller */
} /* namespace ctr */

namespace ctr {
namespace manual {

class AggregateState {
 public:
  AggregateState(nc::net::GraphNodeIndex src, nc::net::GraphNodeIndex dst,
                 const nc::htsim::MatchRuleKey& key)
      : key_(key), src_(src), dst_(dst) {}

  const nc::htsim::MatchRuleKey& key() const { return key_; }

  nc::net::GraphNodeIndex src() const { return src_; }

  nc::net::GraphNodeIndex dst() const { return dst_; }

 private:
  // Key to match at source and endpoints.
  nc::htsim::MatchRuleKey key_;
  nc::net::GraphNodeIndex src_;
  nc::net::GraphNodeIndex dst_;
};

// A container class that constructs a network, per-device TLDR instance and a
// controller. Aggregates and packet generators can be added to the network.
class NetworkContainer {
 public:
  NetworkContainer(std::chrono::milliseconds max_queue_depth,
                   const nc::net::GraphStorage* graph,
                   nc::EventQueue* event_queue);

  // Adds a new aggregate.
  nc::htsim::MatchRuleKey AddAggregate(const AggregateId& id);

  // Adds all devices and links from the graph to the container. Will construct
  // devices using the device factory. If external_internal_observer is not null
  // will add it to all devices.
  void AddElementsFromGraph(
      ctr::controller::DeviceFactory* device_factory,
      nc::htsim::PacketObserver* external_internal_observer = nullptr,
      nc::htsim::PacketObserver* internal_external_observer = nullptr);

  nc::EventQueue* event_queue() { return event_queue_; }

  nc::htsim::Network* network() { return &network_; }

  const nc::net::GraphStorage* graph() { return graph_; }

  std::vector<const nc::htsim::Queue*> queues() const;

  // Installs one or more paths for a given aggregate.
  void InstallPaths(const ctr::AggregateId& id,
                    const std::vector<RouteAndFraction>& routes_and_fractions);

 private:
  nc::htsim::DeviceInterface* AddOrFindDevice(
      const std::string& id, ctr::controller::DeviceFactory* device_factory);

  void SendPathMessages(const AggregateId& id, const ::nc::net::Walk* path);

  nc::htsim::PacketTag TagForPath(const nc::net::Walk* path);

  nc::EventQueue* event_queue_;
  const nc::net::GraphStorage* graph_;

  std::map<AggregateId, AggregateState> id_to_aggregate_state_;
  nc::htsim::Network network_;

  // Max queue size.
  std::chrono::milliseconds max_queue_depth_;

  // Stores devices, indexed by string not by NodeIndex.
  std::map<std::string, std::unique_ptr<nc::htsim::DeviceInterface>> devices_;
  std::vector<std::unique_ptr<nc::htsim::Pipe>> pipes_;

  // Stores queues in the network.
  std::vector<std::unique_ptr<nc::htsim::Queue>> queues_;

  // Boring traffic dies here.
  nc::htsim::DummyPacketHandler dummy_handler_;

  // Each aggregate will be assigned a unique tag, that external packets
  // arriving at the source will be tagged with.
  std::map<AggregateId, nc::htsim::PacketTag> aggregate_id_to_tag;

  // Paths installed.
  std::map<AggregateId, std::vector<const nc::net::Walk*>>
      id_to_paths_installed_;

  std::map<const nc::net::Walk*, nc::htsim::PacketTag,
           nc::net::WalkPtrComparator> path_to_tag_;
};

}  // namespace manual
}  // namespace ctr
#endif
