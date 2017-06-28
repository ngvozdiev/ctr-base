#include "net_mock.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <memory>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "controller.h"

namespace ctr {

void MockDevice::HandlePacketFromPort(nc::htsim::Port* input_port,
                                      nc::htsim::PacketPtr pkt) {
  nc::Unused(input_port);
  CHECK(pkt->five_tuple().ip_dst() == ip_address());
  HandlePacket(std::move(pkt));
}

void MockDevice::HandleStateUpdate(const nc::htsim::SSCPAddOrUpdate& update) {
  AdvanceTime();

  const nc::htsim::MatchRule& rule = update.rule();
  const nc::htsim::MatchRuleKey& key = rule.key();
  AggregateState& state = nc::FindOrDie(states_, key);

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

    auto it_and_bool = state.paths.emplace(
        std::piecewise_construct, std::forward_as_tuple(tag),
        std::forward_as_tuple(path, action->output_port(), tag));
    PathState& path_state = it_and_bool.first->second;
    path_state.fraction = fraction;
  }

  // There should be no paths removed.
  for (const auto& tag_and_path_state : state.paths) {
    CHECK(nc::ContainsKey(tags_in_actions, tag_and_path_state.first));
  }
}

void MockDevice::AdvanceTime() {
  using namespace std::chrono;

  auto now = event_queue_->CurrentTime();
  microseconds now_micros =
      duration_cast<microseconds>(event_queue_->TimeToNanos(now));
  microseconds last_advance_time_micros = duration_cast<microseconds>(
      event_queue_->TimeToNanos(last_advance_time_));
  if (now_micros == last_advance_time_micros) {
    return;
  }

  for (auto& key_and_state : states_) {
    AggregateState& state = key_and_state.second;
    BinSequence& bin_sequence = state.bin_sequence;
    microseconds bin_size = bin_sequence.bin_size();

    // Need to convert microseconds to bin counts. Will also check to make
    // sure they are proper multiples of eachother.
    size_t now_bin_count = now_micros.count() / bin_size.count();
    CHECK(now_micros.count() % bin_size.count() == 0);

    size_t last_advance_bin_count =
        last_advance_time_micros.count() / bin_size.count();
    CHECK(last_advance_time_micros.count() % bin_size.count() == 0);

    BinSequence to_end = bin_sequence.CutFromStart(now_bin_count);
    BinSequence period_sequence = to_end.Offset(last_advance_bin_count);

    for (auto& tag_and_path_state : state.paths) {
      PathState& path_state = tag_and_path_state.second;
      BinSequence split_sequence =
          period_sequence.PreciseSplitOrDie({path_state.fraction})[0];

      std::vector<PcapDataTraceBin> bins =
          split_sequence.AccumulateBins(bin_size);
      for (const auto& bin : bins) {
        path_state.stats.total_bytes_matched += bin.bytes;
        path_state.stats.total_pkts_matched += bin.packets;
        path_state.total_syns += bin.flows_enter;
      }
    }
  }

  last_advance_time_ = now;
}

void MockDevice::ReplyToRequest(const nc::htsim::SSCPStatsRequest& request,
                                nc::htsim::SSCPStatsReply* reply) {
  AdvanceTime();
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
      auto reply = nc::make_unique<SSCPAck>(
          ip_address_, pkt->five_tuple().ip_src(), event_queue_->CurrentTime(),
          add_or_update_message->tx_id());

      LOG(INFO) << "Will TX ACK " << reply->ToString();
      replies_handler_->HandlePacket(std::move(reply));
    }
  } else if (type == SSCPStatsRequest::kSSCPStatsRequestType) {
    SSCPStatsRequest* stats_request_message =
        static_cast<SSCPStatsRequest*>(pkt.get());

    auto reply = nc::make_unique<SSCPStatsReply>(
        ip_address_, pkt->five_tuple().ip_src(), event_queue_->CurrentTime());
    ReplyToRequest(*stats_request_message, reply.get());

    CHECK(replies_handler_) << "Received stats request, but no output handler";
    replies_handler_->HandlePacket(std::move(reply));
  }
}

}  // namespace ctr
