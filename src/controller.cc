#include "controller.h"

#include <iterator>
#include <limits>
#include <string>
#include <tuple>
#include <type_traits>

#include "ncode_common/src/htsim/tcp.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/substitute.h"
#include "routing_system.h"

namespace ctr {
namespace controller {

constexpr nc::net::IPAddress NetworkContainer::kSinkDeviceAddress;
constexpr nc::net::DevicePortNumber NetworkContainer::kDeviceToSinkPort;

void Controller::InitAggregatesAndNodes(
    const std::map<nc::net::GraphNodeIndex, TLDRNode>& tldrs,
    const std::map<AggregateId, AggregateState>& id_to_aggregate_state) {
  tldrs_ = tldrs;
  id_to_aggregate_state_ = id_to_aggregate_state;

  std::map<AggregateId, AggregateHistory> dummy_requests;
  for (const auto& id_and_aggregate_state : id_to_aggregate_state) {
    const AggregateId& id = id_and_aggregate_state.first;
    const AggregateState& aggregate_state = id_and_aggregate_state.second;
    dummy_requests.emplace(id, aggregate_state.initial_history());
  }
  HandleRequest(dummy_requests, std::numeric_limits<uint64_t>::max(), false);
}

Controller::MessageAndNode Controller::GetUpdate(
    std::unique_ptr<nc::htsim::MatchRule> rule, const TLDRNode* node) {
  auto route_update = nc::make_unique<nc::htsim::SSCPAddOrUpdate>(
      controller_ip_, node->tldr_address, event_queue_->CurrentTime(),
      std::move(rule));

  uint64_t tx_id = ++ack_tx_id_;
  route_update->set_tx_id(tx_id);
  return {std::move(route_update), node};
}

void Controller::PopulatePathMessages(const AggregateId& id,
                                      const ::nc::net::Walk* path,
                                      std::vector<MessageAndNode>* out) {
  using namespace nc::net;
  using namespace nc::htsim;

  std::vector<const nc::net::Walk*>& paths_installed =
      id_to_paths_installed_[id];

  if (std::find(paths_installed.begin(), paths_installed.end(), path) !=
      paths_installed.end()) {
    return;
  }

  LOG(INFO) << "Will install path " << path->ToStringNoPorts(*graph_);
  const Links& links = path->links();
  CHECK(links.size() >= 1);

  // The input port of the first update will be the destination port of the
  // first link.
  nc::net::GraphLinkIndex first_link_index = links.front();
  DevicePortNumber input_port = graph_->GetLink(first_link_index)->dst_port();
  PacketTag tag(TagForPathOrDie(path));
  for (size_t i = 1; i < links.size(); ++i) {
    nc::net::GraphLinkIndex link_index = links[i];
    const nc::net::GraphLink* link = graph_->GetLink(link_index);
    const TLDRNode& tldr_node = nc::FindOrDie(tldrs_, link->src());

    auto new_action =
        nc::make_unique<MatchRuleAction>(link->src_port(), tag, 1);

    nc::net::DevicePortNumber input_port_to_match = input_port;
    MatchRuleKey key(tag, input_port_to_match, {FiveTuple::kDefaultTuple});

    auto new_rule = nc::make_unique<nc::htsim::MatchRule>(key);
    new_rule->AddAction(std::move(new_action));
    out->emplace_back(GetUpdate(std::move(new_rule), &tldr_node));
    input_port = link->dst_port();
  }

  // When traffic arrives at the destination it is up to the destination to
  // handle it according to its own routing. Will not install a route at the
  // destination.
  paths_installed.emplace_back(path);
}

void Controller::HandleAck(uint64_t tx_id) {
  std::vector<uint64_t> to_erase;

  // Have to figure out which round this tx belongs to.
  bool found = false;
  for (auto& round_id_and_round_state : round_states_) {
    uint64_t round_id = round_id_and_round_state.first;
    RoundState& round_state = round_id_and_round_state.second;

    if (nc::ContainsKey(round_state.outstanding_tx_ids, tx_id)) {
      found = true;
      round_state.outstanding_tx_ids.erase(tx_id);
      if (!round_state.outstanding_tx_ids.empty()) {
        continue;
      }

      LOG(INFO) << "Collected ACKs from all updates, will commit round "
                << round_id;
      CommitPendingOutput(&round_state);
      to_erase.emplace_back(round_id);
    }
  }

  CHECK(found) << "ACK with unknown tx id";
  CHECK(to_erase.size() < 2) << "Tx id for multiple rounds";
  for (uint64_t round_id : to_erase) {
    round_states_.erase(round_id);
  }
}

void Controller::HandleRequest(
    const std::map<AggregateId, AggregateHistory>& requests, uint64_t round_id,
    bool quick) {
  RoundState& round_state =
      round_states_.emplace(std::piecewise_construct,
                            std::forward_as_tuple(round_id),
                            std::forward_as_tuple(round_id, quick))
          .first->second;
  if (round_state.quick != quick) {
    LOG(ERROR) << "Quick / non-quick mismatch";
    // There is a race between sending regular non-quick requests and sending
    // quick requests, a round can have different quick/non-quick flags seen. In
    // this case we always update it to non-quick.
    round_state.quick = false;
  }

  // Get the requests for this round, optimization will be triggered when
  // requests are collected from all aggregates in a round.
  std::map<AggregateId, AggregateHistory>& histories = round_state.histories;

  for (const auto& id_and_request : requests) {
    const AggregateId& id = id_and_request.first;
    CHECK(nc::ContainsKey(id_to_aggregate_state_, id)) << "Unknown aggregate";

    const AggregateHistory& request = id_and_request.second;
    nc::InsertOrDieNoPrint(&histories, id, request);
  }

  LOG(INFO) << nc::Substitute(
      "Collected requests from $0 / $1 aggregates round $2", histories.size(),
      id_to_aggregate_state_.size(), round_id);

  if (histories.size() == id_to_aggregate_state_.size()) {
    // We heard from all aggregates. Time to re-optimize.
    ReOptimize(&round_state);
  }
}

std::map<AggregateId, AggregateUpdateState> Controller::RoutingToUpdateState(
    const RoutingConfiguration& routing_config) {
  std::map<AggregateId, AggregateUpdateState> out;
  for (const auto& id_and_routes : routing_config.routes()) {
    const AggregateId& id = id_and_routes.first;
    const std::vector<RouteAndFraction>& routes = id_and_routes.second;
    const AggregateState& aggregate_state =
        nc::FindOrDieNoPrint(id_to_aggregate_state_, id);

    AggregateUpdateState update_state(id, aggregate_state.key());
    for (const auto& route_and_fraction : routes) {
      // The output port is the first port in the path.
      const nc::net::Walk* path = route_and_fraction.first;
      const nc::net::Links& links = path->links();
      nc::net::GraphLinkIndex first_link_index = links.front();
      const nc::net::GraphLink* first_link = graph_->GetLink(first_link_index);

      nc::net::DevicePortNumber src_port = first_link->src_port();
      nc::htsim::PacketTag tag(TagForPath(path));
      PathUpdateState path_update_state(tag, src_port, route_and_fraction);
      update_state.paths.emplace_back(path_update_state);
    }

    // TODO: populate the aggregate's limit.
    out.emplace(id, update_state);
  }

  return out;
}

uint32_t Controller::TagForPath(const nc::net::Walk* path) {
  uint32_t* current_tag = nc::FindOrNull(path_to_tag_, path->links());
  if (current_tag != nullptr) {
    return *current_tag;
  }

  uint32_t new_tag = path_to_tag_.size() + 1;
  path_to_tag_.emplace(path->links(), new_tag);
  return new_tag;
}

uint32_t Controller::TagForPathOrDie(const nc::net::Walk* path) const {
  return nc::FindOrDieNoPrint(path_to_tag_, path->links());
}

void Controller::ReOptimize(RoundState* round) {
  std::map<AggregateId, AggregateUpdateState>& pending_output =
      round->pending_output;
  std::set<uint64_t>& outstanding_tx_ids = round->outstanding_tx_ids;

  CHECK(!pending_output.empty());
  CHECK(outstanding_tx_ids.empty());
  nc::EventQueueTime delta = event_queue_->CurrentTime() - last_optimize_time_;
  if (last_optimize_time_ != nc::EventQueueTime::ZeroTime() &&
      delta < event_queue_->ToTime(std::chrono::milliseconds(0))) {
    return;
  }

  LOG(INFO) << "Will re-optimize the network";
  pending_output =
      RoutingToUpdateState(*routing_system_->Update(round->histories));
  LOG(INFO) << "Done optimizing";

  last_optimize_time_ = event_queue_->CurrentTime();

  // We will install paths. InstallPath may modify pending_output_ if all
  // handlers are directly connected and there is no delay. In this case the
  // requests will produce acks which will be handled *before* InstallPath
  // returns. To avoid memory problems we first collect all paths to be
  // installed and then install them.
  std::vector<std::pair<AggregateId, const nc::net::Walk*>> paths_to_install;
  for (const auto& id_and_update_state : pending_output) {
    const AggregateId& id = id_and_update_state.first;
    const AggregateUpdateState& update_state = id_and_update_state.second;

    for (const auto& path_update_state : update_state.paths) {
      const nc::net::Walk* path = path_update_state.route_and_fraction.first;
      paths_to_install.emplace_back(id, path);
    }
  }

  std::vector<MessageAndNode> messages_to_send;
  for (const auto& cookie_and_path : paths_to_install) {
    PopulatePathMessages(cookie_and_path.first, cookie_and_path.second,
                         &messages_to_send);
  }

  // Need to do this in 2 passes - the message may return as soon as it is sent,
  // which will cause problems.
  for (const auto& message_and_node : messages_to_send) {
    outstanding_tx_ids.emplace(message_and_node.message->tx_id());
  }

  for (auto& message_and_node : messages_to_send) {
    LOG(INFO) << "Tx " << message_and_node.message->ToString();
    PacketHandler* to_tldr = message_and_node.node->to_tldr;
    to_tldr->HandlePacket(std::move(message_and_node.message));
  }

  if (outstanding_tx_ids.empty()) {
    // No new paths had to be installed in the network. We can just go ahead
    // and commit the pending output. The pending output may be gone already --
    // the receipt of an ACK during InstallPath (if there are no delays and the
    // handlers are directly connected) may have triggered it already.
    if (!pending_output.empty()) {
      CommitPendingOutput(round);
      round_states_.erase(round->round_id);
    }
  } else {
    // Will wait for the paths to be installed.
    LOG(INFO) << "Need to wait for " << outstanding_tx_ids.size() << " ACKs";
  }
}

void Controller::HandlePacket(::nc::htsim::PacketPtr pkt) {
  using nc::htsim::SSCPAddOrUpdate;
  using nc::htsim::SSCPAck;

  if (pkt->size_bytes() == 0) {
    uint8_t type = pkt->five_tuple().ip_proto().Raw();

    if (type == SSCPAck::kSSCPAckType) {
      SSCPAck* ack_message = static_cast<SSCPAck*>(pkt.get());
      LOG(INFO) << "Rx " << ack_message->ToString();
      HandleAck(ack_message->tx_id());
    } else if (type == TLDRRequest::kTLDRRequestType) {
      TLDRRequest* request_message = static_cast<TLDRRequest*>(pkt.get());
      LOG(INFO) << "Rx " << request_message->ToString();
      HandleRequest(request_message->aggregates(), request_message->round_id(),
                    request_message->quick());
    } else if (type == TLDRTriggerReoptimize::kTLDRTriggerReoptimizeType) {
      for (const auto& src_and_tldr : tldrs_) {
        const TLDRNode* node = &src_and_tldr.second;
        auto message_to_send = ::nc::make_unique<TLDRForceRequest>(
            node->tldr_address, controller_ip_, event_queue_->CurrentTime());
        node->to_tldr->HandlePacket(std::move(message_to_send));
      }
    }

    return;
  }
}

void Controller::CommitPendingOutput(RoundState* round_state) {
  std::map<AggregateId, AggregateUpdateState>& pending_output =
      round_state->pending_output;
  CHECK(!pending_output.empty());

  std::map<const TLDRNode*, std::map<AggregateId, AggregateUpdateState>>
      tldr_node_to_update_map;

  for (auto& id_and_aggregate_state : id_to_aggregate_state_) {
    const AggregateId& id = id_and_aggregate_state.first;
    const AggregateState& aggregate_state = id_and_aggregate_state.second;
    const AggregateUpdateState& update =
        nc::FindOrDieNoPrint(pending_output, id);

    const TLDRNode& node = nc::FindOrDie(tldrs_, aggregate_state.src());
    tldr_node_to_update_map[&node].emplace(id, update);
  }

  for (const auto& tldr_node_ptr_and_update_map : tldr_node_to_update_map) {
    const TLDRNode* node = tldr_node_ptr_and_update_map.first;
    const std::map<AggregateId, AggregateUpdateState>& update_map =
        tldr_node_ptr_and_update_map.second;

    auto message_to_send =
        ::nc::make_unique<TLDRUpdate>(node->tldr_address, controller_ip_,
                                      event_queue_->CurrentTime(), update_map);
    node->to_tldr->HandlePacket(std::move(message_to_send));
  }
}

NetworkContainer::NetworkContainer(const NetworkContainerConfig& config,
                                   const TLDRConfig& tldr_config,
                                   std::unique_ptr<Controller> controller,
                                   const nc::net::GraphStorage* graph,
                                   nc::EventQueue* event_queue)
    : config_(config),
      tldr_config_(tldr_config),
      event_queue_(event_queue),
      graph_(graph),
      controller_(std::move(controller)),
      network_(event_queue->ToTime(std::chrono::milliseconds(10)),
               event_queue) {
  CHECK(tldr_config.ip_src == nc::htsim::kWildIPAddress)
      << "TLDR config should have default source address";
  CHECK(tldr_config.ip_switch_dst == nc::htsim::kWildIPAddress)
      << "TLDR config should have default device address";

  CHECK(config_.min_delay_controller_tldr <= config.max_delay_controller_tldr);
  CHECK(config_.min_delay_tldr_device <= config.max_delay_tldr_device);

  FromGraph();
}

void NetworkContainer::AddAggregate(const AggregateId& id,
                                    AggregateState aggregate_state) {
  id_to_aggregate_state_.emplace(id, aggregate_state);
}

void NetworkContainer::AddGenerator(
    const std::string& input_device, nc::htsim::PacketTag tag,
    std::vector<std::unique_ptr<nc::htsim::BulkPacketSource>> sources) {
  std::string gen_id =
      nc::StrCat("gen_", std::to_string(packet_generators_.size()));

  auto& device = nc::FindOrDie(devices_, input_device);
  nc::htsim::Port* enter_port =
      device->FindOrCreatePort(config_.default_enter_port);

  auto new_gen = nc::make_unique<nc::htsim::BulkPacketGenerator>(
      gen_id, std::move(sources), enter_port, event_queue_);

  // All incoming packets will be tagged with the same tag.
  new_gen->set_default_tag(tag);
  packet_generators_.emplace_back(std::move(new_gen));
}

void NetworkContainer::Init() {
  controller_->InitAggregatesAndNodes(device_id_to_tldr_node_,
                                      id_to_aggregate_state_);
}

nc::htsim::Device* NetworkContainer::AddOrFindDevice(
    const std::string& id,
    std::uniform_int_distribution<size_t>* tldr_device_delay_dist,
    std::uniform_int_distribution<size_t>* controller_tldr_device_dist,
    std::default_random_engine* engine) {
  auto it = devices_.find(id);
  if (it != devices_.end()) {
    return it->second.get();
  }

  std::chrono::milliseconds delay_tldr_device(
      (*tldr_device_delay_dist)(*engine));
  std::chrono::milliseconds delay_controller_tldr(
      (*controller_tldr_device_dist)(*engine));

  size_t n = devices_.size();
  nc::net::IPAddress address =
      nc::net::IPAddress(config_.device_ip_address_base.Raw() + n);
  nc::net::IPAddress tldr_address =
      nc::net::IPAddress(config_.tldr_ip_address_base.Raw() + n);

  auto new_device =
      nc::make_unique<nc::htsim::Device>(id, address, event_queue_);
  network_.AddDevice(new_device.get());

  // Port that accepts traffic from the TLDR.
  nc::htsim::Port* from_tldr =
      new_device->FindOrCreatePort(config_.default_tldr_input_port);

  std::string tldr_id = nc::StrCat(id, "_tldr");
  nc::htsim::PacketHandler* tldr_to_device = from_tldr;
  nc::htsim::PacketHandler* tldr_to_controller = controller_.get();
  if (delay_tldr_device.count()) {
    auto new_pipe = nc::make_unique<nc::htsim::Pipe>(
        tldr_id, id, event_queue_->ToTime(delay_tldr_device), event_queue_);
    new_pipe->Connect(from_tldr);
    tldr_to_device = new_pipe.get();
    pipes_.emplace_back(std::move(new_pipe));
  }

  if (delay_controller_tldr.count()) {
    auto new_pipe = nc::make_unique<nc::htsim::Pipe>(
        tldr_id, "controller", event_queue_->ToTime(delay_controller_tldr),
        event_queue_);
    new_pipe->Connect(controller_.get());
    tldr_to_controller = new_pipe.get();
    pipes_.emplace_back(std::move(new_pipe));
  }

  TLDRConfig tldr_config = tldr_config_;
  tldr_config.ip_src = tldr_address;
  tldr_config.ip_switch_dst = address;

  auto new_tldr =
      nc::make_unique<TLDR>(tldr_id, tldr_config, tldr_to_device,
                            tldr_to_controller, graph_, event_queue_);
  nc::htsim::PacketHandler* device_to_tldr = new_tldr.get();
  nc::htsim::PacketHandler* controller_to_tldr = new_tldr.get();

  if (delay_controller_tldr.count()) {
    auto new_pipe = nc::make_unique<nc::htsim::Pipe>(
        "controller", tldr_id, event_queue_->ToTime(delay_controller_tldr),
        event_queue_);
    new_pipe->Connect(new_tldr.get());
    controller_to_tldr = new_pipe.get();
    pipes_.emplace_back(std::move(new_pipe));
  }

  if (delay_tldr_device.count()) {
    auto new_pipe = nc::make_unique<nc::htsim::Pipe>(
        id, tldr_id, event_queue_->ToTime(delay_tldr_device), event_queue_);
    new_pipe->Connect(new_tldr.get());
    device_to_tldr = new_pipe.get();
    pipes_.emplace_back(std::move(new_pipe));
  }

  device_id_to_tldr_node_.emplace(
      std::piecewise_construct,
      std::forward_as_tuple(graph_->NodeFromStringOrDie(id)),
      std::forward_as_tuple(tldr_address, controller_to_tldr));

  nc::htsim::Device* raw_ptr = new_device.get();
  new_device->set_tx_replies_handler(device_to_tldr);
  new_device->EnableSampling(device_to_tldr, 100);
  devices_.emplace(id, std::move(new_device));
  tldrs_.emplace_back(std::move(new_tldr));

  return raw_ptr;
}

static void AddSrcOrDstBasedRoute(nc::htsim::Device* device,
                                  nc::net::IPAddress address, bool src,
                                  nc::net::DevicePortNumber out_port,
                                  nc::htsim::PacketTag tag) {
  using namespace nc::htsim;

  nc::net::FiveTuple tuple_to_match;
  if (src) {
    tuple_to_match =
        nc::net::FiveTuple(address, kWildIPAddress, kWildIPProto,
                           kWildAccessLayerPort, kWildAccessLayerPort);
  } else {
    tuple_to_match =
        nc::net::FiveTuple(kWildIPAddress, address, kWildIPProto,
                           kWildAccessLayerPort, kWildAccessLayerPort);
  }

  MatchRuleKey key(kWildPacketTag, kWildDevicePortNumber, {tuple_to_match});
  auto action = nc::make_unique<MatchRuleAction>(out_port, tag, 100);
  auto rule = nc::make_unique<MatchRule>(key);
  rule->AddAction(std::move(action));

  auto message =
      nc::make_unique<SSCPAddOrUpdate>(kWildIPAddress, device->ip_address(),
                                       nc::EventQueueTime(0), std::move(rule));
  device->HandlePacket(std::move(message));
}

static void AddDstBasedRoute(nc::htsim::Device* device, nc::net::IPAddress dst,
                             nc::net::DevicePortNumber out_port,
                             nc::htsim::PacketTag tag) {
  AddSrcOrDstBasedRoute(device, dst, false, out_port, tag);
}

static void AddSrcBasedRoute(nc::htsim::Device* device, nc::net::IPAddress src,
                             nc::net::DevicePortNumber out_port,
                             nc::htsim::PacketTag tag) {
  AddSrcOrDstBasedRoute(device, src, true, out_port, tag);
}

static void AddDefaultRoute(nc::htsim::Device* device,
                            nc::net::DevicePortNumber out_port,
                            nc::htsim::PacketTag tag) {
  AddDstBasedRoute(device, nc::htsim::kWildIPAddress, out_port, tag);
}

void NetworkContainer::ConnectToSink(const std::string& device,
                                     nc::net::IPAddress src,
                                     nc::htsim::PacketTag reverse_tag) {
  auto it = devices_.find(device);
  CHECK(it != devices_.end());
  nc::htsim::Device* device_ptr = it->second.get();
  nc::htsim::Device* sink_ptr = GetSinkDevice();

  nc::htsim::Port* port_on_device =
      device_ptr->FindOrCreatePort(config_.default_enter_port);
  nc::htsim::Port* port_on_sink = sink_ptr->NextAvailablePort();
  port_on_sink->Connect(port_on_device);

  AddSrcBasedRoute(device_ptr, src, kDeviceToSinkPort,
                   nc::htsim::kNullPacketTag);
  AddDstBasedRoute(sink_ptr, src, port_on_sink->number(), reverse_tag);
}

std::pair<nc::htsim::Connection*, nc::htsim::Queue*>
NetworkContainer::AddTCPSource(nc::net::IPAddress ip_source,
                               nc::htsim::PacketTag forward_tag,
                               const std::string& source,
                               std::chrono::milliseconds delay,
                               size_t max_queue_size_bytes,
                               nc::net::Bandwidth forward_queue_rate,
                               bool random_queue, double seed) {
  CHECK(!source.empty());

  // Let's first add the new device where the connection will originate.
  std::string id = nc::StrCat("D_", std::to_string(devices_.size()));
  auto new_device =
      nc::make_unique<nc::htsim::Device>(id, ip_source, event_queue_);
  network_.AddDevice(new_device.get());
  AddDefaultRoute(new_device.get(), nc::net::DevicePortNumber(1), forward_tag);

  // Add a unidirectional link between the new device and the source. The source
  // will RX packets via the default_enter_port.
  std::unique_ptr<nc::htsim::Queue> forward_queue;
  if (random_queue) {
    forward_queue = nc::make_unique<nc::htsim::RandomQueue>(
        id, source, forward_queue_rate, max_queue_size_bytes,
        max_queue_size_bytes * 0.5, seed, event_queue_);
    LOG(INFO) << "Added TCP source with RED queue " << max_queue_size_bytes;
  } else {
    forward_queue = nc::make_unique<nc::htsim::FIFOQueue>(
        id, source, forward_queue_rate, max_queue_size_bytes, event_queue_);
    LOG(INFO) << "Added TCP source with FIFO queue " << max_queue_size_bytes;
  }

  auto forward_pipe = nc::make_unique<nc::htsim::Pipe>(
      id, source, event_queue_->ToTime(delay), event_queue_);
  network_.AddLink(forward_queue.get(), forward_pipe.get());

  // Also have to pick a new port on the source and connect it to the new device
  // so that the return path works.
  auto it = devices_.find(source);
  CHECK(it != devices_.end());
  nc::htsim::Device* source_device = it->second.get();
  nc::htsim::Port* port = source_device->NextAvailablePort();
  port->Connect(new_device->FindOrCreatePort(nc::net::DevicePortNumber(1)));
  AddDstBasedRoute(source_device, ip_source, port->number(),
                   nc::htsim::kDefaultTag);

  nc::htsim::Connection* new_connection = new_device->AddTCPGenerator(
      GetSinkDevice()->ip_address(), nc::net::AccessLayerPort(1), kDefaultMSS,
      kDefaultMaxCWND);
  nc::htsim::Queue* forward_queue_raw_ptr = forward_queue.get();

  devices_.emplace(id, std::move(new_device));
  queues_.emplace_back(std::move(forward_queue));
  pipes_.emplace_back(std::move(forward_pipe));
  return {new_connection, forward_queue_raw_ptr};
}

void NetworkContainer::AddDefaultRouteToDummyHandler(const std::string device) {
  auto it = devices_.find(device);
  CHECK(it != devices_.end());
  nc::htsim::Device* device_ptr = it->second.get();

  nc::htsim::Port* port = device_ptr->NextAvailablePort();
  port->Connect(&dummy_handler_);
  AddDefaultRoute(device_ptr, port->number(), nc::htsim::kDefaultTag);
}

void NetworkContainer::FromGraph() {
  std::map<std::string, nc::htsim::Device*> devices;

  std::default_random_engine generator(config_.seed);
  std::uniform_int_distribution<size_t> tldr_device_dist(
      config_.min_delay_tldr_device.count(),
      config_.max_delay_tldr_device.count());
  std::uniform_int_distribution<size_t> controller_tldr_dist(
      config_.min_delay_controller_tldr.count(),
      config_.max_delay_controller_tldr.count());

  size_t links_seed = 0;
  for (nc::net::GraphLinkIndex link : graph_->AllLinks()) {
    const nc::net::GraphLink* link_ptr = graph_->GetLink(link);

    AddOrFindDevice(link_ptr->src_id(), &tldr_device_dist,
                    &controller_tldr_dist, &generator);
    AddOrFindDevice(link_ptr->dst_id(), &tldr_device_dist,
                    &controller_tldr_dist, &generator);

    uint64_t speed_bps = link_ptr->bandwidth().bps();
    uint64_t speed_Bps = speed_bps / 8.0;
    double bytes_per_millisecond = speed_Bps / 1000.0;

    auto new_pipe = nc::make_unique<nc::htsim::Pipe>(*link_ptr, event_queue_);
    uint64_t size = bytes_per_millisecond * config_.max_queue_depth.count();

    std::unique_ptr<nc::htsim::Queue> new_queue;
    if (config_.random_queues) {
      new_queue = nc::make_unique<nc::htsim::RandomQueue>(
          *link_ptr, size, size * 0.5, ++links_seed, event_queue_);
    } else {
      new_queue =
          nc::make_unique<nc::htsim::FIFOQueue>(*link_ptr, size, event_queue_);
    }

    // AddLink will take care of connecting the src/dst ports for us.
    network_.AddLink(new_queue.get(), new_pipe.get(), true);
    queues_.emplace_back(std::move(new_queue));
    pipes_.emplace_back(std::move(new_pipe));
  }
}

nc::htsim::Device* NetworkContainer::GetSinkDevice() {
  if (sink_device_) {
    return sink_device_.get();
  }

  sink_device_ = nc::make_unique<nc::htsim::Device>("Sink", kSinkDeviceAddress,
                                                    event_queue_);

  for (const auto& id_and_device : devices_) {
    nc::htsim::Device* device = id_and_device.second.get();
    nc::htsim::Port* port = device->FindOrCreatePort(kDeviceToSinkPort);
    nc::htsim::Port* sink_port =
        sink_device_->FindOrCreatePort(kDeviceToSinkPort);
    port->Connect(sink_port);
  }

  network_.AddDevice(sink_device_.get());

  return sink_device_.get();
}

}  // namespace controller
}  // namespace ctr
