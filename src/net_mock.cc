#include "net_mock.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "metrics/metrics.h"

static auto* per_link_residuals =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<double, std::string>(
        "link_residuals", "Remaining capacity for a link", "link");

namespace ctr {

void MockDevice::HandlePacketFromPort(nc::htsim::Port* input_port,
                                      nc::htsim::PacketPtr pkt) {
  nc::Unused(input_port);
  CHECK(pkt->five_tuple().ip_dst() == ip_address());
  HandlePacket(std::move(pkt));
}

void MockDevice::HandleStateUpdate(const nc::htsim::SSCPAddOrUpdate& update) {
  mock_network_->AdvanceTime();

  const nc::htsim::MatchRule& rule = update.rule();
  const nc::htsim::MatchRuleKey& key = rule.key();
  AggregateState* state = nc::FindOrNull(states_, key);
  if (state == nullptr) {
    return;
  }

  double total_weight = 0;
  for (const nc::htsim::MatchRuleAction* action : rule.actions()) {
    total_weight += action->weight();
  }

  // Each action in the rule is a separate path.
  std::set<nc::htsim::PacketTag> tags_in_actions;
  for (const nc::htsim::MatchRuleAction* action : rule.actions()) {
    nc::htsim::PacketTag tag = action->tag();
    CHECK(tag.IsNotZero());
    tags_in_actions.emplace(tag);

    const nc::net::Walk* path = controller_->PathForTagOrDie(tag.Raw());
    double fraction = action->weight() / total_weight;

    auto it_and_bool = state->paths.emplace(
        std::piecewise_construct, std::forward_as_tuple(tag),
        std::forward_as_tuple(path, action->output_port(), tag));
    PathState& path_state = it_and_bool.first->second;
    path_state.fraction = fraction;
  }

  // There should be no paths removed.
  for (const auto& tag_and_path_state : state->paths) {
    CHECK(nc::ContainsKey(tags_in_actions, tag_and_path_state.first));
  }
}

void MockNetwork::AdvanceTime() {
  using namespace std::chrono;

  auto now = event_queue_->CurrentTime();
  microseconds now_micros =
      duration_cast<microseconds>(event_queue_->TimeToNanos(now));
  microseconds last_advance_time_micros = duration_cast<microseconds>(
      event_queue_->TimeToNanos(last_advance_time_));
  if (now_micros == last_advance_time_micros) {
    return;
  }

  // Will assume that all bin sizes are the same.
  CHECK(!devices_.empty() && !devices_[0]->states_.empty());
  const MockDevice::AggregateState& aggregate_state =
      devices_[0]->states_.begin()->second;
  microseconds bin_size = aggregate_state.bin_sequence.bin_size();

  // Need to convert microseconds to bin counts. Will also check to make
  // sure they are proper multiples of eachother.
  size_t now_bin_count = now_micros.count() / bin_size.count();
  CHECK(now_micros.count() % bin_size.count() == 0);

  size_t last_advance_bin_count =
      last_advance_time_micros.count() / bin_size.count();
  CHECK(last_advance_time_micros.count() % bin_size.count() == 0);

  // For each link, a BinSequence of all bins that cross the link.
  nc::net::GraphLinkMap<std::unique_ptr<BinSequence>> bins_per_link;

  for (MockDevice* device : devices_) {
    for (auto& key_and_state : device->states_) {
      MockDevice::AggregateState& state = key_and_state.second;
      BinSequence& bin_sequence = state.bin_sequence;

      BinSequence to_end = bin_sequence.CutFromStart(now_bin_count);
      BinSequence period_sequence = to_end.Offset(last_advance_bin_count);

      for (auto& tag_and_path_state : state.paths) {
        MockDevice::PathState& path_state = tag_and_path_state.second;
        BinSequence split_sequence =
            period_sequence.PreciseSplitOrDie({path_state.fraction})[0];
        for (nc::net::GraphLinkIndex link : path_state.path->links()) {
          std::unique_ptr<BinSequence>& bins_for_link = bins_per_link[link];
          if (bins_for_link) {
            bins_for_link->Combine(split_sequence);
          } else {
            bins_for_link = nc::make_unique<BinSequence>(split_sequence);
          }
        }

        std::vector<PcapDataTraceBin> bins =
            split_sequence.AccumulateBins(bin_size);
        for (const auto& bin : bins) {
          path_state.stats.total_bytes_matched += bin.bytes;
          path_state.stats.total_pkts_matched += bin.packets;
          path_state.total_syns += bin.flows_enter;
        }
      }
    }
  }

  for (const auto& link_and_bins : bins_per_link) {
    nc::net::GraphLinkIndex link = link_and_bins.first;
    const nc::net::GraphLink* link_ptr = graph_->GetLink(link);

    BinSequence& bin_sequence = *(*link_and_bins.second);
    nc::net::Bandwidth link_bandwidth = link_ptr->bandwidth();

    std::vector<double> residuals = bin_sequence.Residuals(link_bandwidth);
    for (size_t i = 0; i < residuals.size(); ++i) {
      uint64_t timestamp =
          last_advance_time_.Raw() + i * event_queue_->ToTime(bin_size).Raw();
      per_link_residuals->GetHandle(link_ptr->ToStringNoPorts())
          ->AddValueWithTimestamp(timestamp, residuals[i]);
    }
  }

  last_advance_time_ = now;
}

void MockDevice::ReplyToRequest(const nc::htsim::SSCPStatsRequest& request,
                                nc::htsim::SSCPStatsReply* reply) {
  mock_network_->AdvanceTime();
  for (auto& key_and_state : states_) {
    AggregateState& state = key_and_state.second;

    std::vector<nc::htsim::ActionStats> stats;
    for (auto& tag_and_path : state.paths) {
      PathState& path_state = tag_and_path.second;

      nc::htsim::ActionStats to_add = path_state.stats;
      if (request.include_flow_counts()) {
        to_add.flow_count = path_state.total_syns;
        path_state.total_syns = 0;
      }

      stats.emplace_back(to_add);
    }

    reply->AddStats(key_and_state.first, stats);
  }
}

void MockDevice::HandlePacket(nc::htsim::PacketPtr pkt) {
  using namespace nc::htsim;
  CHECK(pkt->size_bytes() == 0);
  uint8_t type = pkt->five_tuple().ip_proto().Raw();
  if (type == SSCPAddOrUpdate::kSSCPAddOrUpdateType) {
    SSCPAddOrUpdate* add_or_update_message =
        static_cast<SSCPAddOrUpdate*>(pkt.get());
    HandleStateUpdate(*add_or_update_message);
    if (add_or_update_message->tx_id() != SSCPMessage::kNoTxId &&
        replies_handler_ != nullptr) {
      auto reply = nc::GetFreeList<SSCPAck>().New(
          ip_address_, pkt->five_tuple().ip_src(), event_queue_->CurrentTime(),
          add_or_update_message->tx_id());

      LOG(INFO) << "Will TX ACK " << reply->ToString();
      replies_handler_->HandlePacket(std::move(reply));
    }
  } else if (type == SSCPStatsRequest::kSSCPStatsRequestType) {
    SSCPStatsRequest* stats_request_message =
        static_cast<SSCPStatsRequest*>(pkt.get());

    auto reply = nc::GetFreeList<SSCPStatsReply>().New(
        ip_address_, pkt->five_tuple().ip_src(), event_queue_->CurrentTime());
    ReplyToRequest(*stats_request_message, reply.get());

    CHECK(replies_handler_) << "Received stats request, but no output handler";
    replies_handler_->HandlePacket(std::move(reply));
  }
}

void MockSimDevice::HandleStateUpdate(
    const nc::htsim::SSCPAddOrUpdate& update) {
  nc::htsim::MatchRule* rule = update.MutableRule();
  const nc::htsim::MatchRuleKey& key = rule->key();
  AggregateState* state = nc::FindOrNull(states_, key);
  if (state == nullptr) {
    return;
  }

  state->rule = rule;
  double total_weight = 0;
  for (const nc::htsim::MatchRuleAction* action : rule->actions()) {
    total_weight += action->weight();
  }

  // Each action in the rule is a separate path.
  std::set<nc::htsim::PacketTag> tags_in_actions;
  for (const nc::htsim::MatchRuleAction* action : rule->actions()) {
    nc::htsim::PacketTag tag = action->tag();
    CHECK(tag.IsNotZero());
    tags_in_actions.emplace(tag);

    double fraction = action->weight() / total_weight;
    PathState& path_state = state->paths[tag];
    path_state.fraction = fraction;
  }

  // There should be no paths removed.
  for (const auto& tag_and_path_state : state->paths) {
    CHECK(nc::ContainsKey(tags_in_actions, tag_and_path_state.first));
  }
}

static nc::htsim::ActionStats* FindStatsForTagOrDie(
    std::vector<nc::htsim::ActionStats>* stats, nc::htsim::PacketTag tag) {
  for (auto& action_stats : *stats) {
    if (action_stats.tag == tag) {
      return &action_stats;
    }
  }

  LOG(FATAL) << "Could not find action stats";
  return nullptr;
}

void MockSimDevice::PostProcessStats(const nc::htsim::SSCPStatsRequest& request,
                                     nc::htsim::SSCPStatsReply* reply) {
  for (auto& key_and_state : states_) {
    const nc::htsim::MatchRuleKey& key = key_and_state.first;
    AggregateState& state = key_and_state.second;
    std::vector<nc::htsim::ActionStats>& stats_in_reply =
        nc::FindOrDie(reply->stats_mutable(), key);

    for (auto& tag_and_path : state.paths) {
      PathState& path_state = tag_and_path.second;
      if (request.include_flow_counts()) {
        nc::htsim::ActionStats* action_stats =
            FindStatsForTagOrDie(&stats_in_reply, tag_and_path.first);
        action_stats->flow_count = path_state.total_syns;
        path_state.total_syns = 0;
      }
    }
  }
}

void MockSimDevice::HandlePacket(nc::htsim::PacketPtr pkt) {
  using namespace nc::htsim;
  if (pkt->size_bytes() == 0) {
    // The packet is an SSCP message.
    uint8_t type = pkt->five_tuple().ip_proto().Raw();
    if (type == SSCPAddOrUpdate::kSSCPAddOrUpdateType) {
      SSCPAddOrUpdate* add_or_update_message =
          static_cast<SSCPAddOrUpdate*>(pkt.get());
      HandleStateUpdate(*add_or_update_message);
    }
  }

  Device::HandlePacket(std::move(pkt));
}

std::chrono::microseconds MockSimNetwork::GetBinSize() {
  // Will assume that all bin sizes are the same.
  CHECK(!devices_.empty() && !devices_[0]->states_.empty());
  const MockSimDevice::AggregateState& aggregate_state =
      devices_[0]->states_.begin()->second;
  return aggregate_state.bin_sequence.bin_size();
}

void MockSimNetwork::AdvanceTimeToNextBin() {
  using namespace std::chrono;
  microseconds bin_size = GetBinSize();

  for (MockSimDevice* device : devices_) {
    for (auto& key_and_state : device->states_) {
      MockSimDevice::AggregateState& aggregate_state = key_and_state.second;
      BinSequence& bin_sequence = aggregate_state.bin_sequence;

      BinSequence to_end = bin_sequence.CutFromStart(last_bin_count_ + 1);
      BinSequence period_sequence = to_end.Offset(last_bin_count_);

      size_t i = -1;
      for (auto& tag_and_path_state : aggregate_state.paths) {
        MockSimDevice::PathState& path_state = tag_and_path_state.second;
        BinSequence split_sequence =
            period_sequence.PreciseSplitOrDie({path_state.fraction})[0];

        std::vector<PcapDataTraceBin> bins =
            split_sequence.AccumulateBins(bin_size);
        CHECK(bins.size() == 1);
        const PcapDataTraceBin& bin = bins[0];
        if (bin.bytes == 0) {
          continue;
        }

        path_state.total_syns += bin.flows_enter;
        nc::htsim::PacketPtr pkt = GetDummyPacket(bin.bytes);
        const nc::htsim::MatchRuleAction* action =
            aggregate_state.rule->ExplicitChooseOrDie(*pkt, ++i);

        nc::htsim::Port* port = device->FindOrCreatePort(default_enter_port_);
        device->HandlePacketWithAction(port, std::move(pkt), action);
      }
    }
  }

  ++last_bin_count_;
}

nc::htsim::PacketPtr MockSimNetwork::GetDummyPacket(uint32_t size) {
  nc::net::FiveTuple five_tuple(
      nc::net::IPAddress(9919), nc::net::IPAddress(9929), nc::net::kProtoUDP,
      nc::net::AccessLayerPort(9919), nc::net::AccessLayerPort(9929));
  return nc::GetFreeList<nc::htsim::UDPPacket>().New(
      five_tuple, size, event_queue_->CurrentTime());
}

}  // namespace ctr
