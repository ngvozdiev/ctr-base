#include "net_mock.h"

#include <gflags/gflags.h>
#include <ncode/free_list.h>
#include <ncode/logging.h>
#include <ncode/map_util.h>
#include <limits>
#include <numeric>
#include <set>
#include <tuple>

#include "common.h"
#include "metrics/metrics.h"
#include "routing_system.h"

namespace ctr {

DEFINE_bool(precise_splits, false, "If true all splits will be precise");

static auto* link_utilization_metric =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "link_utilization", "Records per-link utilization", "Link source",
            "Link destination");

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
  sub_sequences = state->initial_bin_sequence->PreciseSplitOrDie(fractions);
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
          ++i;
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

NetMock::NetMock(std::map<AggregateId, BinSequence>&& initial_sequences,
                 std::chrono::milliseconds period_duration,
                 std::chrono::milliseconds history_bin_size,
                 RoutingSystem* routing_system)
    : history_bin_size_(history_bin_size),
      initial_sequences_(std::move(initial_sequences)),
      routing_system_(routing_system),
      graph_(routing_system_->graph()) {
  CHECK(!initial_sequences_.empty());
  size_t min_bin_count = std::numeric_limits<size_t>::max();
  std::chrono::milliseconds bin_size = std::chrono::milliseconds::zero();
  for (const auto& id_and_bin_sequence : initial_sequences_) {
    const BinSequence& bin_sequence = id_and_bin_sequence.second;
    if (bin_size == std::chrono::milliseconds::zero()) {
      bin_size = std::chrono::duration_cast<std::chrono::milliseconds>(
          bin_sequence.bin_size());
    } else {
      CHECK(bin_size == bin_sequence.bin_size());
    }

    size_t count = bin_sequence.bin_count();
    min_bin_count = std::min(min_bin_count, count);
  }

  period_duration_bins_ = period_duration.count() / bin_size.count();
  CHECK(period_duration_bins_ > 0);
  period_count_ = min_bin_count / period_duration_bins_;
}

// Generates the input to the system.
std::map<AggregateId, AggregateHistory> NetMock::GenerateInput(
    const std::map<AggregateId, BinSequence>& period_sequences,
    PcapDataBinCache* cache) const {
  std::map<AggregateId, AggregateHistory> input;
  for (const auto& aggregate_and_bins : period_sequences) {
    const AggregateId& aggregate = aggregate_and_bins.first;
    const BinSequence& bins = aggregate_and_bins.second;

    input.emplace(std::piecewise_construct, std::forward_as_tuple(aggregate),
                  std::forward_as_tuple(
                      bins.GenerateHistory(history_bin_size_, 1000, cache)));
  }

  return input;
}

nc::net::GraphLinkMap<std::vector<double>> NetMock::CheckOutput(
    const std::map<AggregateId, BinSequence>& period_sequences,
    const RoutingConfiguration& configuration, PcapDataBinCache* cache) const {
  nc::net::GraphLinkMap<std::unique_ptr<BinSequence>> link_to_bins;

  // First need to figure out which paths cross each link. Will also build a
  // map from paths to aggregates and path indices.
  for (const auto& aggregate_and_routes : configuration.routes()) {
    const AggregateId& aggregate = aggregate_and_routes.first;
    const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;

    std::vector<double> fractions;
    for (const auto& route_and_fraction : routes) {
      fractions.emplace_back(route_and_fraction.second);
    }

    // For each of the aggregate's paths, the bins that go on that path.
    std::vector<std::unique_ptr<BinSequence>> aggregate_split =
        nc::FindOrDieNoPrint(period_sequences, aggregate)
            .PreciseSplitOrDie(fractions);

    for (size_t i = 0; i < routes.size(); ++i) {
      const nc::net::Walk* path = routes[i].first;
      for (nc::net::GraphLinkIndex link : path->links()) {
        std::unique_ptr<BinSequence>& bin_sequence_ptr = link_to_bins[link];
        if (!bin_sequence_ptr) {
          bin_sequence_ptr = aggregate_split[i]->Duplicate();
        } else {
          bin_sequence_ptr->Combine(*aggregate_split[i]);
        }
      }
    }
  }

  nc::net::GraphLinkMap<std::vector<double>> out;
  for (const auto& link_and_bins : link_to_bins) {
    nc::net::GraphLinkIndex link = link_and_bins.first;
    nc::net::Bandwidth rate = graph_->GetLink(link)->bandwidth();
    out[link] = (*link_and_bins.second)->Residuals(rate, cache);
  }

  return out;
}

std::map<AggregateId, BinSequence> NetMock::GetNthPeriod(size_t n) const {
  size_t period_start_bin = n * period_duration_bins_;
  size_t period_end_bin = (n + 1) * period_duration_bins_;

  std::map<AggregateId, BinSequence> out;
  for (const auto& aggregate_and_bins : initial_sequences_) {
    const AggregateId& aggregate = aggregate_and_bins.first;
    const BinSequence& bins = aggregate_and_bins.second;

    std::unique_ptr<BinSequence> to_end = bins.CutFromStart(period_end_bin);
    std::unique_ptr<BinSequence> period_sequence =
        to_end->Offset(period_start_bin);
    out.emplace(std::piecewise_construct, std::forward_as_tuple(aggregate),
                std::forward_as_tuple(period_sequence->traces()));
  }

  return out;
}

std::unique_ptr<RoutingConfiguration> NetMock::InitialOutput(
    PcapDataBinCache* cache) const {
  std::map<AggregateId, BinSequence> zero_period = GetNthPeriod(0);
  std::map<AggregateId, AggregateHistory> input =
      GenerateInput(zero_period, cache);
  return routing_system_->Update(input).routing;
}

static size_t CheckSameSize(
    const nc::net::GraphLinkMap<std::vector<double>>& values) {
  size_t i = 0;
  for (const auto& link_and_values : values) {
    const std::vector<double>& v = *link_and_values.second;
    CHECK(v.size() > 0);
    if (i == 0) {
      i = v.size();
    } else {
      CHECK(i == v.size());
    }
  }

  return i;
}

void NetMock::Run(PcapDataBinCache* cache) {
  std::unique_ptr<RoutingConfiguration> output = InitialOutput(cache);
  size_t timestamp = 0;
  for (size_t i = 0; i < period_count_; ++i) {
    LOG(ERROR) << "Period " << i;

    std::map<AggregateId, BinSequence> period_sequences = GetNthPeriod(i);
    nc::net::GraphLinkMap<std::vector<double>> per_link_residuals =
        CheckOutput(period_sequences, *output, cache);
    size_t num_residuals = CheckSameSize(per_link_residuals);

    for (auto link_and_residuals : per_link_residuals) {
      nc::net::GraphLinkIndex link = link_and_residuals.first;
      std::vector<double>& residuals = *link_and_residuals.second;

      const nc::net::GraphLink* link_ptr = graph_->GetLink(link);
      auto* handle = link_utilization_metric->GetHandle(link_ptr->src_id(),
                                                        link_ptr->dst_id());
      size_t t = timestamp;
      for (double v : residuals) {
        handle->AddValueWithTimestamp(t++, v);
      }
    }

    timestamp += num_residuals;
    std::map<AggregateId, AggregateHistory> input =
        GenerateInput(period_sequences, cache);
    output = std::move(routing_system_->Update(input).routing);
  }
}

}  // namespace ctr
