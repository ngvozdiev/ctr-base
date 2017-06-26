#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <stddef.h>
#include <stdint.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <utility>
#include <vector>

namespace nc {
namespace htsim {
class BulkPacketGenerator;
class BulkPacketSource;
} /* namespace htsim */
} /* namespace ncode */

namespace ctr {
namespace controller {

// A simple component that records per-path packets per period and bits per
// period. It observes all packets that enter the network. It assumes that all
// packets that it sees are tagged with a path's tag. It maintains separate
// metrics for packets with preferential dropping and packets without. Also
// updates the flow recorder.
class InputPacketObserver : public nc::htsim::PacketObserver {
 public:
  InputPacketObserver(nc::net::GraphStorage* path_storage);

  void ObservePacket(const nc::htsim::Packet& pkt);

 private:
  struct PerPathState {
    PerPathState() : bytes_seen(0), pkts_seen(0) {}

    uint64_t bytes_seen;
    uint64_t pkts_seen;
  };

  // Indexed by tag value, and a bool indicating whether or not the path is
  // preferential drop.
  std::map<std::pair<nc::htsim::PacketTag, bool>, PerPathState>
      tag_to_per_path_state_;

  // Needed to relate tags to path strings.
  nc::net::GraphStorage* graph_;
};

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
  // Radius in time from the waypoint, outside this radius the waypoint has no
  // effect, inside the radius values are scaled linearly the closer they are to
  // the waypoint.
  static constexpr std::chrono::milliseconds kWaypointRadius =
      std::chrono::minutes(2);

  AggregateState(nc::net::GraphNodeIndex src, nc::net::GraphNodeIndex dst,
                 const nc::htsim::MatchRuleKey& key,
                 const AggregateHistory& initial_history)
      : key_(key), src_(src), dst_(dst), initial_history_(initial_history) {}

  // This aggregate's effective demand at a given point in time when the
  // waypoint is taken into account.
  nc::net::Bandwidth GetEffectiveDemand(nc::net::Bandwidth demand,
                                        nc::EventQueueTime at,
                                        nc::EventQueue* event_queue) const;

  // If 'now' is past the waypoint (2 * radius + (waypoint_end -
  // waypoint_start)) will move the waypoint to 'now'. If 'now' is within the
  // radius of the waypoint will extend it (move waypoint_end) to now + radius.
  void AddWaypoint(nc::net::Bandwidth target_demand,
                   nc::EventQueue* event_queue);

  const nc::htsim::MatchRuleKey& key() const { return key_; }

  nc::net::GraphNodeIndex src() const { return src_; }

  nc::net::GraphNodeIndex dst() const { return dst_; }

  const AggregateHistory& initial_history() const { return initial_history_; }

 private:
  using Waypoint = std::pair<nc::net::Bandwidth, nc::EventQueueTime>;

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

  Controller(nc::net::IPAddress ip_address, nc::EventQueue* event_queue)
      : controller_ip_(ip_address),
        event_queue_(event_queue),
        ack_tx_id_(0),
        last_optimize_time_(nc::EventQueueTime::ZeroTime()) {}

  // Initializes the state needed for controller operation and bootstraps all
  // aggregates with their initial requests.
  void InitAggregatesAndNodes(
      const std::map<nc::net::GraphNodeIndex, TLDRNode>& tldrs,
      const std::map<uint64_t, AggregateState>& cookie_to_aggregate_state);

  void HandlePacket(::nc::htsim::PacketPtr pkt) override;

  const std::map<uint64_t, AggregateState>& cookie_to_aggregate_state() const {
    return cookie_to_aggregate_state_;
  }

 private:
  struct MessageAndNode {
    std::unique_ptr<nc::htsim::SSCPAddOrUpdate> message;
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
    std::unique_ptr<RoutingConfiguration> pending_output;

    // Per-aggregate histories.
    std::map<uint64_t, AggregateHistory> histories;

    // Path updates that still need to be acknowledged.
    std::set<uint64_t> outstanding_tx_ids;
  };

  // Generates an input suitable for optimization.
  fubar::Input GenerateInput(
      const std::map<uint64_t, tldr::AggregateHistory>& histories,
      const std::map<uint64_t, size_t>& flow_counts,
      const std::map<uint64_t, nc::net::Bandwidth>& rates) const;

  // Updates the target of aggregates that are oversubscribed to improve
  // fairness.
  void AddWaypointsToOversubscribedAggregates(
      const fubar::Output& output,
      const std::map<uint64_t, tldr::AggregateHistory>& histories);

  MessageAndNode GetUpdate(std::unique_ptr<nc::htsim::MatchRule> rule,
                           const TLDRNode* node);

  // Populates a vector with messages to install the given path if the path is
  // not installed already.
  void PopulatePathMessages(uint64_t cookie, const ::nc::net::GraphPath* path,
                            std::vector<MessageAndNode>* out);

  void HandleAck(uint64_t tx_id);

  void HandleRequest(const std::map<uint64_t, tldr::AggregateHistory>& requests,
                     uint64_t round_id, bool quick);

  // Re-optimizes the network based on each aggregate's history.
  void ReOptimize(RoundState* round_state);

  // Optimizes the network with progressively larger aggregate scale so that
  // queuing is below a threshold.
  std::unique_ptr<fubar::Output> IterativeOptimize(
      const std::map<uint64_t, tldr::AggregateHistory>& histories, bool quick,
      uint64_t round_id);

  void CommitPendingOutput(RoundState* round_state);

  fubar::Output FlipFilter(const fubar::Output& output);

  // All traffic originating from the controller will have this ip.
  const nc::net::IPAddress controller_ip_;

  // A map from a node id in the graph to the channel that can be used to send
  // messages to the TLDR instance that controls it.
  std::map<nc::net::GraphNodeIndex, TLDRNode> tldrs_;

  // Stores per-aggregate information.
  std::map<uint64_t, AggregateState> cookie_to_aggregate_state_;

  // Keeps track of the paths currently installed in the network. This set is
  // updated as soon as messages are sent to install the paths, not when they
  // are acknowledged. Indexed by aggregate cookie.
  std::map<uint64_t, std::vector<const nc::net::GraphPath*>>
      cookie_to_paths_installed_;

  // The event queue.
  nc::EventQueue* event_queue_;

  // Id for ACKs.
  uint64_t ack_tx_id_;

  // Records the time the optimizer last ran. Useful when avoiding running the
  // optimizer too often.
  nc::EventQueueTime last_optimize_time_;

  // Records the most recent optimization output.
  std::unique_ptr<fubar::Output> last_optimize_output_;

  // Per-optimization round states.
  std::map<uint64_t, RoundState> round_states_;
};

struct NetworkContainerConfig {
  NetworkContainerConfig(nc::net::IPAddress device_ip_address_base,
                         nc::net::IPAddress tldr_ip_address_base,
                         nc::net::DevicePortNumber default_enter_port,
                         nc::net::DevicePortNumber default_tldr_input_port,
                         std::chrono::milliseconds max_queue_depth,
                         bool random_queues)
      : device_ip_address_base(device_ip_address_base),
        tldr_ip_address_base(tldr_ip_address_base),
        default_enter_port(default_enter_port),
        default_tldr_input_port(default_tldr_input_port),
        min_delay_tldr_device(0),
        max_delay_tldr_device(0),
        min_delay_controller_tldr(0),
        max_delay_controller_tldr(0),
        max_queue_depth(max_queue_depth),
        seed(0),
        random_queues(random_queues) {}

  // Devices in the network will have addresses assigned starting with this
  // address.
  nc::net::IPAddress device_ip_address_base;

  // TLDR instances in the network will have addresses assigned starting with
  // this address.
  nc::net::IPAddress tldr_ip_address_base;

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
};

// A container class that constructs a network, per-device TLDR instance and a
// controller. Aggregates and packet generators can be added to the network.
class NetworkContainer {
 public:
  static constexpr size_t kDefaultMSS = 1500;
  static constexpr size_t kDefaultMaxCWND = 2000000;
  static constexpr nc::net::IPAddress kSinkDeviceAddress =
      nc::net::IPAddress(1);
  static constexpr nc::net::DevicePortNumber kDeviceToSinkPort =
      nc::net::DevicePortNumber(9995);

  NetworkContainer(const nc::net::PBNet& net_pb,
                   const NetworkContainerConfig& config,
                   const tldr::TLDRConfig& tldr_config,
                   std::unique_ptr<Controller> controller,
                   nc::net::GraphStorage* storage, nc::EventQueue* event_queue);

  // Adds a new aggregate.
  void AddAggregate(uint64_t cookie, AggregateState aggregate_state);

  // Adds a number of sources attached at 'input_device'. All packets will be
  // tagged with the given tag.
  void AddGenerator(
      const std::string& input_device, nc::htsim::PacketTag tag,
      std::vector<std::unique_ptr<nc::htsim::BulkPacketSource>> sources);

  // Connects a port on a device directly to the universal sink. Also installs a
  // route at the sink to direct traffic from the sink to a given destination
  // back to the device.
  void ConnectToSink(const std::string& device, nc::net::IPAddress src,
                     nc::htsim::PacketTag reverse_tag);

  // Adds a default route from the given device to a dummy handler. Useful if
  // you worry about dropping packets at the device.
  void AddDefaultRouteToDummyHandler(const std::string device);

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
      const std::string& source, std::chrono::milliseconds delay,
      size_t max_queue_size_bytes, uint64_t forward_queue_rate_bps,
      bool random_queue, double seed, bool important = true);

  void Init();

  nc::EventQueue* event_queue() { return event_queue_; }

  const NetworkContainerConfig& config() const { return config_; }
  const tldr::TLDRConfig& tldr_config() const { return tldr_config_; }

  nc::htsim::Network* network() { return &network_; }

  const nc::net::GraphStorage* storage() { return storage_; }

 private:
  nc::htsim::Device* AddOrFindDevice(
      const std::string& id,
      std::uniform_int_distribution<size_t>* tldr_device_delay_dist,
      std::uniform_int_distribution<size_t>* controller_tldr_device_dist,
      std::default_random_engine* engine);

  nc::htsim::Device* GetSinkDevice();

  void FromPB(const nc::net::PBNet& net_pb);

  // Configuration of the container.
  const NetworkContainerConfig config_;

  // Base TLDR config. All TLDR instances will be constructed with this config,
  // but with different addresses.
  const tldr::TLDRConfig tldr_config_;

  nc::EventQueue* event_queue_;
  nc::net::GraphStorage* storage_;

  // The controller config and the controller.
  std::map<uint64_t, AggregateState> cookie_to_aggregate_state_;
  std::unique_ptr<Controller> controller_;

  nc::htsim::Network network_;
  std::map<nc::net::GraphNodeIndex, TLDRNode> device_id_to_tldr_node_;
  std::map<std::string, std::unique_ptr<nc::htsim::Device>> devices_;
  std::vector<std::unique_ptr<nc::htsim::Pipe>> pipes_;
  std::vector<std::unique_ptr<nc::htsim::Queue>> queues_;

  // The packet sources.
  std::vector<std::unique_ptr<nc::htsim::BulkPacketGenerator>>
      packet_generators_;

  // TLDR instances, one for each device.
  std::vector<std::unique_ptr<tldr::TLDR>> tldrs_;

  // TCP connections are terminated here.
  std::unique_ptr<nc::htsim::Device> sink_device_;

  // Boring traffic dies here.
  nc::htsim::DummyPacketHandler dummy_handler_;

  // Input observer that sees packets as they enter the network.
  InputPacketObserver input_packet_observer_;
};

// A controller that always uses the shortest path.
class SPController : public Controller {
 public:
  using Controller::Controller;

  std::unique_ptr<fubar::Output> Optimize(
      const fubar::Input& input, const fubar::PathMap* path_map) override;
};

// A controller that runs FUBAR. All aggregates will have the same cost.
class FUBARController : public Controller {
 public:
  using Controller::Controller;

  std::unique_ptr<fubar::Output> Optimize(
      const fubar::Input& input, const fubar::PathMap* path_map) override;
};

class MinMaxController : public Controller {
 public:
  using Controller::Controller;

  std::unique_ptr<fubar::Output> Optimize(
      const fubar::Input& input, const fubar::PathMap* path_map) override;
};

class B4Controller : public Controller {
 public:
  using Controller::Controller;

  std::unique_ptr<fubar::Output> Optimize(
      const fubar::Input& input, const fubar::PathMap* path_map) override;
};

}  // namespace controller
}  // namespace ctr
#endif
