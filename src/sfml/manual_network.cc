#include "manual_network.h"

#include <stddef.h>
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <tuple>

#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/free_list.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "../net_mock.h"

namespace ctr {
namespace manual {

constexpr nc::net::DevicePortNumber NetworkContainer::kDefaultEnterPort;
static constexpr nc::net::DevicePortNumber kDefaultDeviceBaseAddress(2000);

NetworkContainer::NetworkContainer(std::chrono::milliseconds max_queue_depth,
                                   const nc::net::GraphStorage* graph,
                                   nc::EventQueue* event_queue)
    : event_queue_(event_queue),
      graph_(graph),
      network_(event_queue->ToTime(std::chrono::hours(1)), event_queue),
      max_queue_depth_(max_queue_depth) {}

std::vector<const nc::htsim::Queue*> NetworkContainer::queues() const {
  std::vector<const nc::htsim::Queue*> queues_raw;
  for (const auto& src_dst_and_queue_ptr : queues_) {
    queues_raw.emplace_back(src_dst_and_queue_ptr.second.get());
  }

  return queues_raw;
}

// Returns a flow key that matches all traffic entering a port.
static nc::htsim::MatchRuleKey GetKey(nc::net::DevicePortNumber input_port,
                                      nc::htsim::PacketTag tag) {
  using namespace nc::htsim;
  nc::net::FiveTuple five_tuple_to_match(kWildIPAddress, kWildIPAddress,
                                         kWildIPProto, kWildAccessLayerPort,
                                         kWildAccessLayerPort);
  MatchRuleKey key(tag, input_port, {five_tuple_to_match});
  return key;
}

nc::htsim::MatchRuleKey NetworkContainer::AddAggregate(
    const AggregateId& aggregate_id) {
  AggregateState* current_state =
      nc::FindOrNull(id_to_aggregate_state_, aggregate_id);
  if (current_state != nullptr) {
    return current_state->key();
  }

  // Need a new tag for the aggregate.
  nc::htsim::PacketTag aggregate_tag(aggregate_id_to_tag.size() + 1);
  aggregate_id_to_tag.emplace(aggregate_id, aggregate_tag);
  nc::htsim::MatchRuleKey aggregate_key =
      GetKey(kDefaultEnterPort, aggregate_tag);

  id_to_aggregate_state_.emplace(
      std::piecewise_construct, std::forward_as_tuple(aggregate_id),
      std::forward_as_tuple(aggregate_id.src(), aggregate_id.dst(),
                            aggregate_key));
  return aggregate_key;
}

nc::htsim::DeviceInterface* NetworkContainer::AddOrFindDevice(
    const std::string& id, ctr::controller::DeviceFactory* device_factory) {
  auto it = devices_.find(id);
  if (it != devices_.end()) {
    return it->second.get();
  }

  size_t n = devices_.size();
  nc::net::IPAddress address =
      nc::net::IPAddress(kDefaultDeviceBaseAddress.Raw() + n);

  std::unique_ptr<nc::htsim::DeviceInterface> new_device =
      device_factory->NewDevice(id, address, event_queue_);
  network_.AddDevice(new_device.get());

  nc::htsim::DeviceInterface* raw_ptr = new_device.get();
  devices_.emplace(id, std::move(new_device));

  return raw_ptr;
}

void NetworkContainer::AddElementsFromGraph(
    ctr::controller::DeviceFactory* device_factory,
    nc::htsim::PacketObserver* external_internal_observer,
    nc::htsim::PacketObserver* internal_external_observer) {
  for (nc::net::GraphLinkIndex link : graph_->AllLinks()) {
    const nc::net::GraphLink* link_ptr = graph_->GetLink(link);

    nc::htsim::DeviceInterface* src_device =
        AddOrFindDevice(link_ptr->src_id(), device_factory);
    nc::htsim::DeviceInterface* dst_device =
        AddOrFindDevice(link_ptr->dst_id(), device_factory);
    if (external_internal_observer) {
      src_device->AddExternalInternalObserver(external_internal_observer);
      dst_device->AddExternalInternalObserver(external_internal_observer);
    }

    if (internal_external_observer) {
      src_device->AddInternalExternalObserver(internal_external_observer);
      dst_device->AddInternalExternalObserver(internal_external_observer);
    }

    uint64_t speed_bps = link_ptr->bandwidth().bps();
    uint64_t speed_Bps = speed_bps / 8.0;
    double bytes_per_millisecond = speed_Bps / 1000.0;

    auto new_pipe = nc::make_unique<nc::htsim::Pipe>(*link_ptr, event_queue_);
    uint64_t size = bytes_per_millisecond * max_queue_depth_.count();

    std::unique_ptr<nc::htsim::Queue> new_queue =
        nc::make_unique<nc::htsim::FIFOQueue>(*link_ptr, size, event_queue_);

    // AddLink will take care of connecting the src/dst ports for us.
    network_.AddLink(new_queue.get(), new_pipe.get(), link_ptr->src_id(),
                     link_ptr->dst_id(), link_ptr->src_port(),
                     link_ptr->dst_port(), true);
    queues_.emplace(std::make_pair(link_ptr->src_id(), link_ptr->dst_id()),
                    std::move(new_queue));
    pipes_.emplace_back(std::move(new_pipe));
  }
}

nc::htsim::PacketTag NetworkContainer::TagForPath(const nc::net::Walk* path) {
  nc::htsim::PacketTag* current_tag = nc::FindOrNull(path_to_tag_, path);
  if (current_tag != nullptr) {
    return *current_tag;
  }

  nc::htsim::PacketTag new_tag = nc::htsim::PacketTag(path_to_tag_.size() + 1);
  path_to_tag_.emplace(path, new_tag);
  return new_tag;
}

void NetworkContainer::InstallPaths(
    const ctr::AggregateId& id,
    const std::vector<RouteAndFraction>& routes_and_fractions) {
  for (const auto& route_and_fraction : routes_and_fractions) {
    SendPathMessages(id, route_and_fraction.first);
  }

  const nc::htsim::MatchRuleKey& key =
      nc::FindOrDieNoPrint(id_to_aggregate_state_, id).key();
  auto match_rule = nc::make_unique<nc::htsim::MatchRule>(key);
  CHECK(!routes_and_fractions.empty());
  for (const RouteAndFraction& route_and_fraction : routes_and_fractions) {
    double fraction = route_and_fraction.second;
    const nc::net::Walk* path = route_and_fraction.first;

    size_t weight = 1000 * fraction;
    nc::htsim::PacketTag tag = nc::FindOrDieNoPrint(path_to_tag_, path);

    const nc::net::GraphLink* first_link =
        graph_->GetLink(path->links().front());
    nc::net::DevicePortNumber src_port = first_link->src_port();

    auto action =
        nc::make_unique<nc::htsim::MatchRuleAction>(src_port, tag, weight);
    match_rule->AddAction(std::move(action));
  }

  const std::string& src_id = graph_->GetNode(id.src())->id();
  nc::htsim::DeviceInterface* device_iface =
      nc::FindOrDie(devices_, src_id).get();
  auto route_update = nc::GetFreeList<nc::htsim::SSCPAddOrUpdate>().New(
      nc::net::IPAddress(9999), device_iface->ip_address(),
      event_queue_->CurrentTime(), std::move(match_rule));
  device_iface->HandlePacket(std::move(route_update));
}

void NetworkContainer::SendPathMessages(const AggregateId& id,
                                        const ::nc::net::Walk* path) {
  using namespace nc::net;
  using namespace nc::htsim;

  std::vector<const nc::net::Walk*>& paths_installed =
      id_to_paths_installed_[id];

  if (std::find(paths_installed.begin(), paths_installed.end(), path) !=
      paths_installed.end()) {
    LOG(INFO) << "Path already installed " << path->ToStringNoPorts(*graph_);
    return;
  }

  LOG(INFO) << "Will install path " << path->ToStringNoPorts(*graph_);
  const Links& links = path->links();
  CHECK(links.size() >= 1);

  // Will first install state in the tail N-1 hops of the path. The input port
  // of the first update will be the destination port of the first link.
  nc::net::GraphLinkIndex first_link_index = links.front();
  DevicePortNumber input_port = graph_->GetLink(first_link_index)->dst_port();
  PacketTag tag = TagForPath(path);
  for (size_t i = 1; i < links.size(); ++i) {
    nc::net::GraphLinkIndex link_index = links[i];
    const nc::net::GraphLink* link = graph_->GetLink(link_index);
    const std::string& src_id = graph_->GetNode(link->src())->id();
    nc::htsim::DeviceInterface* device_iface =
        nc::FindOrDie(devices_, src_id).get();

    auto new_action =
        nc::make_unique<MatchRuleAction>(link->src_port(), tag, 1);
    nc::net::DevicePortNumber input_port_to_match = input_port;
    MatchRuleKey key(tag, input_port_to_match, {FiveTuple::kDefaultTuple});

    // Will directly update the route with no delay.
    auto new_rule = nc::make_unique<nc::htsim::MatchRule>(key);
    new_rule->AddAction(std::move(new_action));
    auto route_update = nc::GetFreeList<nc::htsim::SSCPAddOrUpdate>().New(
        nc::net::IPAddress(9999), device_iface->ip_address(),
        event_queue_->CurrentTime(), std::move(new_rule));
    device_iface->HandlePacket(std::move(route_update));

    input_port = link->dst_port();
  }
  paths_installed.emplace_back(path);
}

}  // namespace manual
}  // namespace ctr
