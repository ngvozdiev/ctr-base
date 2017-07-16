#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/animator.h"
#include "ncode_common/src/htsim/bulk_gen.h"
#include "ncode_common/src/htsim/flow_driver.h"
#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/htsim/network.h"
#include "ncode_common/src/htsim/packet.h"
#include "ncode_common/src/htsim/queue.h"
#include "ncode_common/src/net/net_common.h"
#include "common.h"
#include "tldr.h"

namespace ctr {
class RoutingSystem;
} /* namespace ctr */

namespace ctr {
namespace controller {

struct TLDRNode {
  TLDRNode(nc::net::IPAddress tldr_address, nc::htsim::PacketHandler* to_tldr)
      : tldr_address(tldr_address), to_tldr(to_tldr) {}

  // Address that messages sent to the TLDR that controls this node will
  // have.
  nc::net::IPAddress tldr_address;

  // Handles messages to TLDR.
  nc::htsim::PacketHandler* to_tldr;
};

class AggregateState {
 public:
  AggregateState(nc::net::GraphNodeIndex src, nc::net::GraphNodeIndex dst,
                 const nc::htsim::MatchRuleKey& key,
                 const AggregateHistory& initial_history)
      : key_(key), src_(src), dst_(dst), initial_history_(initial_history) {}

  const nc::htsim::MatchRuleKey& key() const { return key_; }

  nc::net::GraphNodeIndex src() const { return src_; }

  nc::net::GraphNodeIndex dst() const { return dst_; }

  const AggregateHistory& initial_history() const { return initial_history_; }

  AggregateHistory& initial_history() { return initial_history_; }

 private:
  // Key to match at source and endpoints.
  nc::htsim::MatchRuleKey key_;
  nc::net::GraphNodeIndex src_;
  nc::net::GraphNodeIndex dst_;

  // The initial rate and flow count for this aggregate -- used when
  // bootstrapping the controller.
  AggregateHistory initial_history_;
};

// The controller proceeds in rounds. Each round starts when it has heard recent
// requests from all devices. It then performs optimization and installs paths
// if needed. After paths are installed it replies to the requests with updates,
// which trigger the newly-installed paths to be used.
class Controller : public ::nc::htsim::PacketHandler {
 public:
  using RateAndFlowCount = std::pair<nc::net::Bandwidth, size_t>;

  Controller(nc::net::IPAddress ip_address, RoutingSystem* routing_system,
             nc::EventQueue* event_queue, const nc::net::GraphStorage* graph)
      : controller_ip_(ip_address),
        event_queue_(event_queue),
        ack_tx_id_(0),
        last_optimize_time_(nc::EventQueueTime::ZeroTime()),
        routing_system_(routing_system),
        graph_(graph) {}

  // Initializes the state needed for controller operation and bootstraps all
  // aggregates with their initial requests.
  void InitAggregatesAndNodes(
      const std::map<nc::net::GraphNodeIndex, TLDRNode>& tldrs,
      const std::map<AggregateId, AggregateState>& id_to_aggregate_state);

  void HandlePacket(::nc::htsim::PacketPtr pkt) override;

  const std::map<AggregateId, AggregateState>& id_to_aggregate_state() const {
    return id_to_aggregate_state_;
  }

  const nc::net::Walk* PathForTagOrDie(nc::htsim::PacketTag tag) const;

  const nc::net::Walk* PathForTagOrNull(nc::htsim::PacketTag tag) const;

  const nc::net::GraphStorage* graph() const { return graph_; }

  const nc::EventQueue* event_queue() const { return event_queue_; }

 private:
  struct MessageAndNode {
    nc::FreeList<nc::htsim::SSCPAddOrUpdate>::Pointer message;
    const TLDRNode* node;
  };

  // State for a single optimization round.
  struct RoundState {
    RoundState(uint64_t round_id, bool quick)
        : round_id(round_id), quick(quick) {}

    // The id of the round.
    uint64_t round_id;

    // If this round is quick only paths currently in the network will be used.
    bool quick;

    // The result of an optimization pass is stored here while the controller is
    // waiting for path updates to be acknowledged.
    std::map<AggregateId, AggregateUpdateState> pending_output;

    // Per-aggregate histories.
    std::map<AggregateId, AggregateHistory> histories;

    // Path updates that still need to be acknowledged.
    std::set<uint64_t> outstanding_tx_ids;
  };

  MessageAndNode GetUpdate(std::unique_ptr<nc::htsim::MatchRule> rule,
                           const TLDRNode* node);

  // Populates a vector with messages to install the given path if the path is
  // not installed already.
  void PopulatePathMessages(const AggregateId& id, const ::nc::net::Walk* path,
                            std::vector<MessageAndNode>* out);

  void HandleAck(uint64_t tx_id);

  void HandleRequest(const std::map<AggregateId, AggregateHistory>& requests,
                     uint64_t round_id, bool quick);

  // Re-optimizes the network based on each aggregate's history.
  void ReOptimize(RoundState* round_state);

  void CommitPendingOutput(RoundState* round_state);

  // Generates AggregateUpdateStates from a RoutingConfiguration.
  std::map<AggregateId, AggregateUpdateState> RoutingToUpdateState(
      const RoutingConfiguration& routing_config,
      const CompetingAggregates& competing_aggregates);

  // Generates new tags for paths or returns old ones.
  nc::htsim::PacketTag TagForPath(const nc::net::Walk* path);

  nc::htsim::PacketTag TagForPathOrDie(const nc::net::Walk* path) const;

  // All traffic originating from the controller will have this ip.
  const nc::net::IPAddress controller_ip_;

  // A map from a node id in the graph to the channel that can be used to send
  // messages to the TLDR instance that controls it.
  std::map<nc::net::GraphNodeIndex, TLDRNode> tldrs_;

  // Stores per-aggregate information.
  std::map<AggregateId, AggregateState> id_to_aggregate_state_;

  // Keeps track of the paths currently installed in the network. This set is
  // updated as soon as messages are sent to install the paths, not when they
  // are acknowledged. Indexed by aggregate id.
  std::map<AggregateId, std::vector<const nc::net::Walk*>>
      id_to_paths_installed_;

  // The event queue.
  nc::EventQueue* event_queue_;

  // Id for ACKs.
  uint64_t ack_tx_id_;

  // Records the time the optimizer last ran. Useful when avoiding running the
  // optimizer too often.
  nc::EventQueueTime last_optimize_time_;

  // Per-optimization round states.
  std::map<uint64_t, RoundState> round_states_;

  // The routing system actually does the prediction and optimization.
  RoutingSystem* routing_system_;

  struct WalkPtrComparator {
    bool operator()(const nc::net::Walk* lhs, const nc::net::Walk* rhs) const {
      return *lhs < *rhs;
    }
  };

  // Per-path tags.
  std::map<const nc::net::Walk*, nc::htsim::PacketTag, WalkPtrComparator>
      path_to_tag_;
  std::map<nc::htsim::PacketTag, const nc::net::Walk*> tag_to_path_;

  const nc::net::GraphStorage* graph_;
};

struct NetworkContainerConfig {
  NetworkContainerConfig(nc::net::IPAddress device_ip_address_base,
                         nc::net::IPAddress tldr_ip_address_base,
                         nc::net::IPAddress tcp_group_source_address_base,
                         nc::net::DevicePortNumber default_enter_port,
                         nc::net::DevicePortNumber default_tldr_input_port,
                         std::chrono::milliseconds max_queue_depth,
                         std::chrono::milliseconds tcp_rto_timer_period,
                         bool random_queues)
      : device_ip_address_base(device_ip_address_base),
        tldr_ip_address_base(tldr_ip_address_base),
        tcp_group_source_ip_address_base(tcp_group_source_address_base),
        default_enter_port(default_enter_port),
        default_tldr_input_port(default_tldr_input_port),
        min_delay_tldr_device(0),
        max_delay_tldr_device(0),
        min_delay_controller_tldr(0),
        max_delay_controller_tldr(0),
        max_queue_depth(max_queue_depth),
        seed(0),
        random_queues(random_queues),
        tcp_rto_timer_period(tcp_rto_timer_period) {}

  // Devices in the network will have addresses assigned starting with this
  // address.
  nc::net::IPAddress device_ip_address_base;

  // TLDR instances in the network will have addresses assigned starting with
  // this address.
  nc::net::IPAddress tldr_ip_address_base;

  // The addresses of all devices in TCP groups will have this address as their
  // base.
  nc::net::IPAddress tcp_group_source_ip_address_base;

  // Traffic that enters the network will do so via this port on the first hop.
  nc::net::DevicePortNumber default_enter_port;

  // Port that TLDR instances use to talk to devices.
  nc::net::DevicePortNumber default_tldr_input_port;

  // Min/max delay between the switch and the TLDR instance that controls it.
  std::chrono::milliseconds min_delay_tldr_device;
  std::chrono::milliseconds max_delay_tldr_device;

  // Min/max delay between the controller and the TLDR instance.
  std::chrono::milliseconds min_delay_controller_tldr;
  std::chrono::milliseconds max_delay_controller_tldr;

  // Queues will be at most this deep.
  std::chrono::milliseconds max_queue_depth;

  // Seed for the delay values between controller/switch/TLDR.
  size_t seed;

  // If true will use RED for all queues.
  bool random_queues;

  // How often to run the RTO timer. Set to a very high value to switch off.
  std::chrono::milliseconds tcp_rto_timer_period;

  // Config shared by all TCP sources.
  nc::htsim::TCPSourceConfig tcp_config;
};

class DeviceFactory {
 public:
  virtual ~DeviceFactory(){};

  virtual std::unique_ptr<nc::htsim::DeviceInterface> NewDevice(
      const std::string& id, nc::net::IPAddress address,
      nc::EventQueue* event_queue) = 0;
};

class RateKeyFrame {
 public:
  RateKeyFrame(std::chrono::milliseconds at, size_t flow_count,
               nc::net::Bandwidth rate)
      : at_(at), flow_count_(flow_count), rate_(rate) {}

  // The time for this key frame.
  std::chrono::milliseconds at() const { return at_; }

  // Number of flow affected by the key frame.
  size_t flow_count() const { return flow_count_; }

  // What the rate of those flows should be at time 'at'.
  nc::net::Bandwidth rate() const { return rate_; }

 private:
  std::chrono::milliseconds at_;
  size_t flow_count_;
  nc::net::Bandwidth rate_;
};

// A bunch of TCP flows. Each flow will have a separate access link, each access
// link will have a rate based on a series of key frames, allowing to change the
// access links' rates as the animation progresses. The flows will have
// data added at either uniform or exponential intervals. The amount of data
// added is either uniform or exponential.
class TCPFlowGroup {
 public:
  TCPFlowGroup(size_t flow_count, std::chrono::milliseconds min_external_delay,
               std::chrono::milliseconds max_external_delay,
               size_t mean_object_size_bytes,
               std::chrono::milliseconds mean_wait_time)
      : flow_count_(flow_count),
        min_external_delay_(min_external_delay),
        max_external_delay_(max_external_delay),
        access_link_rate_spread_(0.0),
        random_access_link_queue_(false),
        mean_object_size_bytes_(mean_object_size_bytes),
        mean_wait_time_(mean_wait_time),
        mean_object_size_fixed_(false),
        mean_wait_time_fixed_(false),
        initial_time_offset_(std::chrono::milliseconds::zero()) {}

  // A constant offset to be applied to when data is added to each TCP flow's
  // connection.
  std::chrono::milliseconds initial_time_offset() const {
    return initial_time_offset_;
  }

  void set_initial_time_offset(std::chrono::milliseconds initial_time_offset) {
    initial_time_offset_ = initial_time_offset;
  }

  // If set to true the same amount of bytes (equal to 'mean_object_size_bytes')
  // will be added to all TCP flows' connections. If false the data added will
  // be exponentially distributed around the mean.
  void set_mean_object_size_fixed(bool mean_object_size_fixed) {
    mean_object_size_fixed_ = mean_object_size_fixed;
  }

  // Same as above, but for wait times between when data is added.
  void set_mean_wait_time_fixed(bool mean_wait_time_fixed) {
    mean_wait_time_fixed_ = mean_wait_time_fixed;
  }

  // If true will use RED on the queues.
  void set_random_access_link_queue(bool random_access_link_queue) {
    random_access_link_queue_ = random_access_link_queue;
  }

  // A value between 0 and 1, if 0 each flow's access link will have the value
  // from the key frame. If spread is 1 each flow's access link value will be
  // uniformly picked between 0 and 2 x the value from the keyframe.
  void set_access_link_rate_spread(double rate_spread) {
    access_link_rate_spread_ = rate_spread;
  }

  // Adds a new keyframe. Keyframes dictate the rates of flows' access links. At
  // least one keyframe is needed.
  void AddKeyFrame(std::chrono::milliseconds at, size_t flow_count,
                   nc::net::Bandwidth rate) {
    key_frames_.emplace_back(at, flow_count, rate);
  }

  size_t flow_count() const { return flow_count_; }

  const std::vector<RateKeyFrame>& key_frames() const { return key_frames_; }

  std::chrono::milliseconds max_external_delay() const {
    return max_external_delay_;
  }

  size_t mean_object_size_bytes() const { return mean_object_size_bytes_; }

  bool mean_object_size_fixed() const { return mean_object_size_fixed_; }

  std::chrono::milliseconds mean_wait_time() const { return mean_wait_time_; }

  bool mean_wait_time_fixed() const { return mean_wait_time_fixed_; }

  std::chrono::milliseconds min_external_delay() const {
    return min_external_delay_;
  }

  bool random_access_link_queue() const { return random_access_link_queue_; }

  double access_link_rate_spread() const { return access_link_rate_spread_; }

 private:
  size_t flow_count_;
  std::chrono::milliseconds min_external_delay_;
  std::chrono::milliseconds max_external_delay_;
  std::vector<RateKeyFrame> key_frames_;
  double access_link_rate_spread_;

  bool random_access_link_queue_;
  size_t mean_object_size_bytes_;
  std::chrono::milliseconds mean_wait_time_;
  bool mean_object_size_fixed_;
  bool mean_wait_time_fixed_;

  std::chrono::milliseconds initial_time_offset_;
};

// A container class that constructs a network, per-device TLDR instance and a
// controller. Aggregates and packet generators can be added to the network.
class NetworkContainer {
 public:
  static constexpr nc::net::IPAddress kSinkDeviceAddress =
      nc::net::IPAddress(1);
  static constexpr nc::net::DevicePortNumber kDeviceToSinkPort =
      nc::net::DevicePortNumber(9995);

  NetworkContainer(const NetworkContainerConfig& config,
                   const TLDRConfig& tldr_config,
                   const nc::net::GraphStorage* graph, Controller* controller,
                   nc::EventQueue* event_queue);

  // Adds a new aggregate.
  nc::htsim::MatchRuleKey AddAggregate(const AggregateId& id,
                                       const AggregateHistory& initial_history);

  // Adds a number of sources attached at 'input_device'. All packets will be
  // tagged with the given tag.
  void AddGenerator(
      const std::string& input_device, nc::htsim::PacketTag tag,
      std::vector<std::unique_ptr<nc::htsim::BulkPacketSource>> sources);

  // Adds a number of TCP sources, one for each flow in a flow group. Attached
  // to each source will be a flow driver, which will add data to it according
  // to the distributions in 'flow_group'.
  void AddTCPFlowGroup(const AggregateId& id, const TCPFlowGroup& flow_group);

  // Adds all devices and links from the graph to the container. Will construct
  // devices using the device factory. If external_internal_observer is not null
  // will add it to all devices.
  void AddElementsFromGraph(
      DeviceFactory* device_factory,
      nc::htsim::PacketObserver* external_internal_observer = nullptr);

  // Inializes the controller with per-aggregate state.
  void InitAggregatesInController();

  nc::EventQueue* event_queue() { return event_queue_; }

  const NetworkContainerConfig& config() const { return config_; }
  const TLDRConfig& tldr_config() const { return tldr_config_; }

  nc::htsim::Network* network() { return &network_; }

  const nc::net::GraphStorage* graph() { return graph_; }

  std::vector<const nc::htsim::Queue*> internal_queues() const;

  const std::vector<nc::htsim::TCPSource*>& flow_group_tcp_sources() const {
    return flow_group_tcp_sources_;
  }

  Controller* controller() { return controller_; }

 private:
  // Adds a default route from the given device to a dummy handler. Useful if
  // you worry about dropping packets at the device.
  void AddDefaultRouteToDummyHandler(nc::htsim::DeviceInterface* device_ptr);

  // Connects a port on a device directly to the universal sink. Also installs a
  // route at the sink to direct traffic from the sink to a given destination
  // back to the device.
  void ConnectToSink(nc::net::GraphNodeIndex dst, nc::net::IPAddress src,
                     nc::htsim::PacketTag reverse_tag);

  // Adds a TCP source that will send data packets between a device one hop
  // before 'source' and 'destination'. The TCP stream will originate at its own
  // device that will be connected to the source via a link (queue and pipe)
  // with the given delay. Returns a pair with the connection and the forward
  // queue. The rate of the queue can be animated, and the connection managed by
  // a FlowDriver. The forward queue will be RED with the given max size and
  // seed. The drop threshold will be set to 50% of the queue size. The new
  // device creates will have an ip of 'ip_source'.
  std::pair<nc::htsim::Connection*, nc::htsim::Queue*> AddTCPSource(
      nc::net::IPAddress ip_source, nc::htsim::PacketTag forward_tag,
      nc::net::GraphNodeIndex source, std::chrono::milliseconds delay,
      size_t max_queue_size_bytes, nc::net::Bandwidth forward_queue_rate,
      bool random_queue, double seed);

  nc::htsim::DeviceInterface* AddOrFindDevice(
      const std::string& id,
      std::uniform_int_distribution<size_t>* tldr_device_delay_dist,
      std::uniform_int_distribution<size_t>* controller_tldr_device_dist,
      std::default_random_engine* engine, DeviceFactory* device_factory);

  nc::htsim::DeviceInterface* GetSinkDevice();

  // Configuration of the container.
  const NetworkContainerConfig config_;

  // Base TLDR config. All TLDR instances will be constructed with this config,
  // but with different addresses.
  const TLDRConfig tldr_config_;

  nc::EventQueue* event_queue_;
  const nc::net::GraphStorage* graph_;

  // The controller config and the controller.
  std::map<AggregateId, AggregateState> id_to_aggregate_state_;
  Controller* controller_;

  nc::htsim::Network network_;
  std::map<nc::net::GraphNodeIndex, TLDRNode> device_id_to_tldr_node_;

  // Stores devices, indexed by string not by NodeIndex, because there will be
  // some "fake" devices added to originate TCP connections.
  std::map<std::string, std::unique_ptr<nc::htsim::DeviceInterface>> devices_;
  std::vector<std::unique_ptr<nc::htsim::Pipe>> pipes_;

  // Stores queues in the network. They are separated into internal and external
  // (added to flow groups).
  std::vector<std::unique_ptr<nc::htsim::Queue>> internal_queues_;
  std::vector<std::unique_ptr<nc::htsim::Queue>> external_queues_;

  // The packet sources.
  std::vector<std::unique_ptr<nc::htsim::BulkPacketGenerator>>
      packet_generators_;

  // TLDR instances, one for each device.
  std::vector<std::unique_ptr<TLDR>> tldrs_;

  // TCP connections are terminated here.
  std::unique_ptr<nc::htsim::DeviceInterface> sink_device_;

  // Boring traffic dies here.
  nc::htsim::DummyPacketHandler dummy_handler_;

  std::vector<std::unique_ptr<nc::htsim::FlowDriver>> flow_drivers_;

  nc::htsim::AnimationContainer animation_container_;

  uint64_t seed_gen_;

  // Each aggregate will be assigned a unique tag, that external packets
  // arriving at the source will be tagged with.
  std::map<AggregateId, nc::htsim::PacketTag> aggregate_id_to_tag;

  // Records all TCP sources created as part of flow groups. Owned by the
  // respective devices.
  std::vector<nc::htsim::TCPSource*> flow_group_tcp_sources_;
};

}  // namespace controller
}  // namespace ctr
#endif
