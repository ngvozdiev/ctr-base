#include "gflags/gflags.h"
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

namespace ctr {

DEFINE_bool(precise_splits, false, "If true all splits will be precise");

void MockSimDevice::HandleStateUpdate(
    const nc::htsim::SSCPAddOrUpdate& update) {
  nc::htsim::MatchRule* rule = update.MutableRule();
  std::vector<const nc::htsim::MatchRuleAction*> rule_actions = rule->actions();

  const nc::htsim::MatchRuleKey& key = rule->key();
  AggregateState* state = nc::FindOrNull(states_, key);
  if (state == nullptr) {
    return;
  }

  state->rule = rule;
  double total_weight = 0;
  for (const nc::htsim::MatchRuleAction* action : rule_actions) {
    total_weight += action->weight();
  }

  std::vector<double> fractions;
  for (const nc::htsim::MatchRuleAction* action : rule_actions) {
    double fraction = action->weight() / total_weight;
    fractions.emplace_back(fraction);
  }

  std::vector<std::unique_ptr<BinSequence>> sub_sequences;
  if (FLAGS_precise_splits) {
    sub_sequences = state->initial_bin_sequence->PreciseSplitOrDie(fractions);
  } else {
    sub_sequences = state->initial_bin_sequence->SplitOrDie(fractions);
  }
  CHECK(sub_sequences.size() == fractions.size());

  // Each action in the rule is a separate path.
  std::set<nc::htsim::PacketTag> tags_in_actions;
  for (size_t i = 0; i < rule_actions.size(); ++i) {
    const nc::htsim::MatchRuleAction* action = rule_actions[i];

    nc::htsim::PacketTag tag = action->tag();
    CHECK(tag.IsNotZero());
    tags_in_actions.emplace(tag);

    PathState& path_state = state->paths[tag];
    path_state.bin_sequence = std::move(sub_sequences[i]);
    path_state.bins.clear();
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
        action_stats->flow_count =
            GetFlowCountFromSyns(path_state.syns.GetValues());
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
  for (const MockSimDevice* device : devices_) {
    if (device->states_.empty()) {
      continue;
    }

    const MockSimDevice::AggregateState& aggregate_state =
        device->states_.begin()->second;
    return aggregate_state.initial_bin_sequence->bin_size();
  }

  LOG(FATAL) << "No devices with aggregate state";
  return std::chrono::microseconds(0);
}

void MockSimNetwork::PrefetchBins(MockSimDevice::PathState* path_state) {
  CHECK(path_state->bin_sequence);
  const BinSequence& bin_sequence = *path_state->bin_sequence;
  std::unique_ptr<BinSequence> to_end =
      bin_sequence.CutFromStart(last_bin_count_ + kPrefetchSize);
  std::unique_ptr<BinSequence> period_sequence =
      to_end->Offset(last_bin_count_);
  std::vector<TrimmedPcapDataTraceBin> bins =
      period_sequence->AccumulateBins(GetBinSize(), &bin_cache_);
  CHECK(bins.size() == kPrefetchSize) << bins.size() << " vs " << kPrefetchSize;

  path_state->bins = std::move(bins);
  path_state->bins_cached_from = last_bin_count_;
}

void MockSimNetwork::AdvanceTimeToNextBin() {
  using namespace std::chrono;

  for (MockSimDevice* device : devices_) {
    for (auto& key_and_state : device->states_) {
      MockSimDevice::AggregateState& aggregate_state = key_and_state.second;

      size_t i = -1;
      for (auto& tag_and_path_state : aggregate_state.paths) {
        MockSimDevice::PathState& path_state = tag_and_path_state.second;
        if (!path_state.bin_sequence) {
          continue;
        }

        const std::vector<TrimmedPcapDataTraceBin>& bins = path_state.bins;
        size_t bin_index = last_bin_count_ + 1 - path_state.bins_cached_from;
        if (bin_index >= bins.size()) {
          PrefetchBins(&path_state);
          bin_index = last_bin_count_ + 1 - path_state.bins_cached_from;
          CHECK(bin_index < bins.size());
        }
        const TrimmedPcapDataTraceBin& bin = bins[bin_index];
        if (bin.bytes == 0) {
          continue;
        }

        path_state.syns.AddValue(bin.flows_enter);
        nc::htsim::PacketPtr pkt = GetDummyPacket(bin.bytes);
        const nc::htsim::MatchRuleAction* action =
            aggregate_state.rule->ExplicitChooseOrDie(*pkt, ++i);
        CHECK(action->weight() > 0) << "Action with zero weight chosen "
                                    << action->ToString() << " rule "
                                    << action->parent_rule()->ToString();

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

uint64_t GetFlowCountFromSyns(const std::vector<uint64_t>& syns) {
  uint64_t total = std::accumulate(syns.begin(), syns.end(), 0ul);
  if (total == 0) {
    return 1ul;
  }

  return std::max(static_cast<uint64_t>(1), total / syns.size());
}

}  // namespace ctr
