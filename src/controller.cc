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
  auto route_update = nc::GetFreeList<nc::htsim::SSCPAddOrUpdate>().New(
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
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(routing_config.demands(), id);
    nc::net::Bandwidth demand = demand_and_flow_count.first;

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

    // TODO: This should be higher (if there is slack capacity) or lower (if
    // traffic does not fit).
    update_state.limit = demand;

    out.emplace(id, update_state);
  }

  return out;
}

uint32_t Controller::TagForPath(const nc::net::Walk* path) {
  uint32_t* current_tag = nc::FindOrNull(path_to_tag_, path);
  if (current_tag != nullptr) {
    return *current_tag;
  }

  uint32_t new_tag = path_to_tag_.size() + 1;
  path_to_tag_.emplace(path, new_tag);
  return new_tag;
}

uint32_t Controller::TagForPathOrDie(const nc::net::Walk* path) const {
  return nc::FindOrDieNoPrint(path_to_tag_, path);
}

void Controller::ReOptimize(RoundState* round) {
  std::map<AggregateId, AggregateUpdateState>& pending_output =
      round->pending_output;
  std::set<uint64_t>& outstanding_tx_ids = round->outstanding_tx_ids;

  CHECK(pending_output.empty());
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
    LOG(ERROR) << message_and_node.message->ToString() << " node "
               << message_and_node.node->tldr_address;
    outstanding_tx_ids.emplace(message_and_node.message->tx_id());
  }

  uint64_t round_id = round->round_id;
  for (auto& message_and_node : messages_to_send) {
    LOG(INFO) << "Tx " << message_and_node.message->ToString();
    PacketHandler* to_tldr = message_and_node.node->to_tldr;
    to_tldr->HandlePacket(std::move(message_and_node.message));
  }

  if (!nc::ContainsKey(round_states_, round_id)) {
    // Round was already committed by HandlePacket immediately returning.
    return;
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
        auto message_to_send = nc::GetFreeList<TLDRForceRequest>().New(
            node->tldr_address, controller_ip_, event_queue_->CurrentTime());
        node->to_tldr->HandlePacket(std::move(message_to_send));
      }
    }

    return;
  }
}

const nc::net::Walk* Controller::PathForTagOrDie(uint32_t tag) const {
  for (const auto& path_and_tag : path_to_tag_) {
    const nc::net::Walk* path = path_and_tag.first;
    uint32_t current_tag = path_and_tag.second;
    if (current_tag == tag) {
      return path;
    }
  }

  LOG(FATAL) << "Unable to find tag for path";
  return nullptr;
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

    auto message_to_send = nc::GetFreeList<TLDRUpdate>().New(
        node->tldr_address, controller_ip_, event_queue_->CurrentTime(),
        update_map);
    node->to_tldr->HandlePacket(std::move(message_to_send));
  }
}

NetworkContainer::NetworkContainer(const NetworkContainerConfig& config,
                                   const TLDRConfig& tldr_config,
                                   const nc::net::GraphStorage* graph,
                                   Controller* controller,
                                   nc::EventQueue* event_queue)
    : config_(config),
      tldr_config_(tldr_config),
      event_queue_(event_queue),
      graph_(graph),
      controller_(controller),
      network_(event_queue->ToTime(config.tcp_rto_timer_period), event_queue),
      animation_container_("AnimationContainer", std::chrono::milliseconds(10),
                           event_queue),
      seed_gen_(1.0) {
  CHECK(tldr_config.ip_src == nc::htsim::kWildIPAddress)
      << "TLDR config should have default source address";
  CHECK(tldr_config.ip_switch_dst == nc::htsim::kWildIPAddress)
      << "TLDR config should have default device address";

  CHECK(config_.min_delay_controller_tldr <= config.max_delay_controller_tldr);
  CHECK(config_.min_delay_tldr_device <= config.max_delay_tldr_device);
}

std::vector<const nc::htsim::Queue*> NetworkContainer::internal_queues() const {
  std::vector<const nc::htsim::Queue*> queues_raw;
  for (const auto& queue_ptr : internal_queues_) {
    queues_raw.emplace_back(queue_ptr.get());
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

static void CombineHistories(AggregateHistory* h1, const AggregateHistory& h2) {
  h1->set_flow_count(h1->flow_count() + h2.flow_count());
  size_t bin_count = h1->bins().size();
  if (h2.bins().size() == 1) {
    for (size_t i = 0; i < bin_count; ++i) {
      h1->SetBin(i, h1->bins()[i] + h2.bins()[0]);
    }
    return;
  }

  CHECK(h1->bins().size() == h2.bins().size()) << "History bin count mismatch "
                                               << h1->bins().size() << " vs "
                                               << h2.bins().size();
  CHECK(h1->bin_size() == h2.bin_size()) << "Bin size mismatch "
                                         << h1->bin_size().count() << " vs "
                                         << h2.bin_size().count();
  for (size_t i = 0; i < bin_count; ++i) {
    h1->SetBin(i, h1->bins()[i] + h2.bins()[i]);
  }
}

nc::htsim::MatchRuleKey NetworkContainer::AddAggregate(
    const AggregateId& aggregate_id, const AggregateHistory& initial_history) {
  AggregateState* current_state =
      nc::FindOrNull(id_to_aggregate_state_, aggregate_id);
  if (current_state != nullptr) {
    // Will just combine the histories.
    CombineHistories(&current_state->initial_history(), initial_history);
    return current_state->key();
  }

  // Need a new tag for the aggregate.
  nc::htsim::PacketTag aggregate_tag(aggregate_id_to_tag.size() + 1);
  aggregate_id_to_tag.emplace(aggregate_id, aggregate_tag);
  nc::htsim::MatchRuleKey aggregate_key =
      GetKey(config_.default_enter_port, aggregate_tag);

  id_to_aggregate_state_.emplace(
      std::piecewise_construct, std::forward_as_tuple(aggregate_id),
      std::forward_as_tuple(aggregate_id.src(), aggregate_id.dst(),
                            aggregate_key, initial_history));
  return aggregate_key;
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

void NetworkContainer::InitAggregatesInController() {
  controller_->InitAggregatesAndNodes(device_id_to_tldr_node_,
                                      id_to_aggregate_state_);
}

nc::htsim::DeviceInterface* NetworkContainer::AddOrFindDevice(
    const std::string& id,
    std::uniform_int_distribution<size_t>* tldr_device_delay_dist,
    std::uniform_int_distribution<size_t>* controller_tldr_device_dist,
    std::default_random_engine* engine, DeviceFactory* device_factory) {
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

  std::unique_ptr<nc::htsim::DeviceInterface> new_device =
      device_factory->NewDevice(id, address, event_queue_);
  network_.AddDevice(new_device.get());

  // Port that accepts traffic from the TLDR.
  nc::htsim::Port* from_tldr =
      new_device->FindOrCreatePort(config_.default_tldr_input_port);

  std::string tldr_id = nc::StrCat(id, "_tldr");
  nc::htsim::PacketHandler* tldr_to_device = from_tldr;
  nc::htsim::PacketHandler* tldr_to_controller = controller_;
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
    new_pipe->Connect(controller_);
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

  nc::htsim::DeviceInterface* raw_ptr = new_device.get();
  new_device->set_tx_replies_handler(device_to_tldr);
  devices_.emplace(id, std::move(new_device));
  tldrs_.emplace_back(std::move(new_tldr));

  return raw_ptr;
}

static void AddSrcOrDstBasedRoute(nc::htsim::DeviceInterface* device,
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

  auto message = nc::GetFreeList<SSCPAddOrUpdate>().New(
      kWildIPAddress, device->ip_address(), nc::EventQueueTime(0),
      std::move(rule));
  device->HandlePacket(std::move(message));
}

static void AddDstBasedRoute(nc::htsim::DeviceInterface* device,
                             nc::net::IPAddress dst,
                             nc::net::DevicePortNumber out_port,
                             nc::htsim::PacketTag tag) {
  AddSrcOrDstBasedRoute(device, dst, false, out_port, tag);
}

static void AddSrcBasedRoute(nc::htsim::DeviceInterface* device,
                             nc::net::IPAddress src,
                             nc::net::DevicePortNumber out_port,
                             nc::htsim::PacketTag tag) {
  AddSrcOrDstBasedRoute(device, src, true, out_port, tag);
}

static void AddDefaultRoute(nc::htsim::DeviceInterface* device,
                            nc::net::DevicePortNumber out_port,
                            nc::htsim::PacketTag tag) {
  AddDstBasedRoute(device, nc::htsim::kWildIPAddress, out_port, tag);
}

void NetworkContainer::ConnectToSink(nc::net::GraphNodeIndex device_index,
                                     nc::net::IPAddress src,
                                     nc::htsim::PacketTag reverse_tag) {
  const std::string& device_id = graph_->GetNode(device_index)->id();
  nc::htsim::DeviceInterface& device =
      nc::FindSmartPtrOrDie(devices_, device_id);
  nc::htsim::DeviceInterface* sink_ptr = GetSinkDevice();

  nc::htsim::Port* port_on_device =
      device.FindOrCreatePort(config_.default_enter_port);
  nc::htsim::Port* port_on_sink = sink_ptr->NextAvailablePort();
  port_on_sink->Connect(port_on_device);

  AddSrcBasedRoute(&device, src, kDeviceToSinkPort, nc::htsim::kNullPacketTag);
  AddDstBasedRoute(sink_ptr, src, port_on_sink->number(), reverse_tag);
}

std::pair<nc::htsim::Connection*, nc::htsim::Queue*>
NetworkContainer::AddTCPSource(nc::net::IPAddress ip_source,
                               nc::htsim::PacketTag forward_tag,
                               nc::net::GraphNodeIndex source,
                               std::chrono::milliseconds delay,
                               size_t max_queue_size_bytes,
                               nc::net::Bandwidth forward_queue_rate,
                               bool random_queue, double seed) {
  const nc::net::GraphNode* node_ptr = graph_->GetNode(source);
  const std::string& src_id = node_ptr->id();

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
        id, src_id, forward_queue_rate, max_queue_size_bytes,
        max_queue_size_bytes * 0.5, seed, event_queue_);
    LOG(INFO) << "Added TCP source with RED queue " << max_queue_size_bytes;
  } else {
    forward_queue = nc::make_unique<nc::htsim::FIFOQueue>(
        id, src_id, forward_queue_rate, max_queue_size_bytes, event_queue_);
    LOG(INFO) << "Added TCP source with FIFO queue " << max_queue_size_bytes;
  }

  auto forward_pipe = nc::make_unique<nc::htsim::Pipe>(
      id, src_id, event_queue_->ToTime(delay), event_queue_);
  network_.AddLink(forward_queue.get(), forward_pipe.get(), id, src_id,
                   nc::net::DevicePortNumber(1), config_.default_enter_port);

  // Also have to pick a new port on the source and connect it to the new device
  // so that the return path works.
  nc::htsim::DeviceInterface& source_device =
      nc::FindSmartPtrOrDie(devices_, src_id);
  nc::htsim::Port* port = source_device.NextAvailablePort();
  port->Connect(new_device->FindOrCreatePort(nc::net::DevicePortNumber(1)));
  AddDstBasedRoute(&source_device, ip_source, port->number(),
                   nc::htsim::kDefaultTag);

  nc::htsim::TCPSource* new_connection = new_device->AddTCPGenerator(
      config_.tcp_config, GetSinkDevice()->ip_address(),
      nc::net::AccessLayerPort(1));
  flow_group_tcp_sources_.emplace_back(new_connection);
  nc::htsim::Queue* forward_queue_raw_ptr = forward_queue.get();

  devices_.emplace(id, std::move(new_device));
  external_queues_.emplace_back(std::move(forward_queue));
  pipes_.emplace_back(std::move(forward_pipe));
  return {new_connection, forward_queue_raw_ptr};
}

void NetworkContainer::AddTCPFlowGroup(const AggregateId& id,
                                       const TCPFlowGroup& flow_group) {
  using namespace std::chrono;
  nc::net::GraphNodeIndex src = id.src();
  nc::net::GraphNodeIndex dst = id.dst();

  std::default_random_engine generator(++seed_gen_);
  std::uniform_int_distribution<uint64_t> time_distribution(
      flow_group.min_external_delay().count(),
      flow_group.max_external_delay().count());

  double spread = flow_group.access_link_rate_spread();
  CHECK(spread <= 1 && spread >= 0);

  CHECK(flow_group.key_frames().size() > 0);
  for (size_t i = 0; i < flow_group.flow_count(); ++i) {
    nc::net::IPAddress src_address(
        config_.tcp_group_source_ip_address_base.Raw() + flow_drivers_.size());

    uint64_t max_rate_bps = 0;
    std::vector<nc::htsim::KeyFrame> key_frames;
    for (const RateKeyFrame& rate_key_frame : flow_group.key_frames()) {
      uint64_t mean_rate_bps =
          rate_key_frame.rate().bps() / flow_group.flow_count();

      std::uniform_int_distribution<uint64_t> rate_distribution(
          (1 - spread) * mean_rate_bps, (1 + spread) * mean_rate_bps);
      double rate_bps = rate_distribution(generator);
      if (rate_bps > max_rate_bps) {
        max_rate_bps = rate_bps;
      }

      key_frames.push_back({rate_key_frame.at(), rate_bps});
    }
    CHECK(key_frames.size() > 0);

    nc::net::Bandwidth start_rate = flow_group.key_frames().front().rate();
    milliseconds delay = milliseconds(time_distribution(generator));

    // The AggregateId object has a convenient method to get the SP delay.
    ctr::AggregateId aggregate_id(src, dst);
    milliseconds sp_delay =
        duration_cast<milliseconds>(aggregate_id.GetSPDelay(*graph_));

    // Will size the queue according to the highest rate -- this may cause
    // some flows to have too deep queues if they start off very small.
    double total_delay_sec = (delay.count() + sp_delay.count()) / 1000.0;
    size_t queue_size_bytes = max_rate_bps / 8.0 * total_delay_sec;

    // Will assume that there are already aggregates added for both the forward
    // and the reverse of the connection.
    nc::htsim::PacketTag forward_tag =
        nc::FindOrDieNoPrint(aggregate_id_to_tag, AggregateId(src, dst));
    nc::htsim::PacketTag reverse_tag =
        nc::FindOrDieNoPrint(aggregate_id_to_tag, AggregateId(dst, src));

    nc::htsim::Connection* connection;
    nc::htsim::Queue* queue;
    std::tie(connection, queue) = AddTCPSource(
        src_address, forward_tag, src, delay, queue_size_bytes, start_rate,
        flow_group.random_access_link_queue(), ++seed_gen_);
    ConnectToSink(dst, src_address, reverse_tag);

    auto animator =
        nc::make_unique<nc::htsim::LinearAnimator>(key_frames, false, queue);
    animation_container_.AddAnimator(std::move(animator));

    auto object_sizes_and_wait_times_gen =
        nc::make_unique<nc::htsim::DefaultObjectSizeAndWaitTimeGenerator>(
            flow_group.mean_object_size_bytes(),
            flow_group.mean_object_size_fixed(), flow_group.mean_wait_time(),
            flow_group.mean_wait_time_fixed(), ++seed_gen_, event_queue_);
    object_sizes_and_wait_times_gen->set_constant_delay(
        flow_group.initial_time_offset());

    std::string driver_id = nc::StrCat(connection->id(), "_driver");
    auto flow_driver = nc::make_unique<nc::htsim::FeedbackLoopFlowDriver>(
        driver_id, std::move(object_sizes_and_wait_times_gen), event_queue_);

    flow_driver->ConnectionAttached(connection);
    flow_drivers_.emplace_back(std::move(flow_driver));
  }
}

void NetworkContainer::AddDefaultRouteToDummyHandler(
    nc::htsim::DeviceInterface* device_ptr) {
  nc::htsim::Port* port = device_ptr->NextAvailablePort();
  port->Connect(&dummy_handler_);
  AddDefaultRoute(device_ptr, port->number(), nc::htsim::kDefaultTag);
}

void NetworkContainer::AddElementsFromGraph(DeviceFactory* device_factory) {
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
                    &controller_tldr_dist, &generator, device_factory);
    AddOrFindDevice(link_ptr->dst_id(), &tldr_device_dist,
                    &controller_tldr_dist, &generator, device_factory);

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
    network_.AddLink(new_queue.get(), new_pipe.get(), link_ptr->src_id(),
                     link_ptr->dst_id(), link_ptr->src_port(),
                     link_ptr->dst_port(), true);
    internal_queues_.emplace_back(std::move(new_queue));
    pipes_.emplace_back(std::move(new_pipe));
  }

  for (auto& id_and_device : devices_) {
    nc::htsim::DeviceInterface* device_ptr = id_and_device.second.get();
    AddDefaultRouteToDummyHandler(device_ptr);
  }
}

nc::htsim::DeviceInterface* NetworkContainer::GetSinkDevice() {
  if (sink_device_) {
    return sink_device_.get();
  }

  sink_device_ = nc::make_unique<nc::htsim::Device>("Sink", kSinkDeviceAddress,
                                                    event_queue_);

  for (const auto& id_and_device : devices_) {
    nc::htsim::DeviceInterface* device = id_and_device.second.get();
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
