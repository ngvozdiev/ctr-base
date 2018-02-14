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

DEFINE_bool(mean_hints, false,
            "If true each aggregate mean level will be set to the ideal one "
            "for the next history.");
DEFINE_bool(
    triggered_optimization, true,
    "If true will convolve aggregates to trigger optimization if needed at the "
    "end of each period.");

static auto* link_utilization_metric =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "link_utilization", "Records per-link utilization (bytes per bin)",
            "Link source", "Link destination");

static auto* queue_size_metric =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "queue_size", "Records per-link queue size (in bytes)",
            "Link source", "Link destination");

static auto* link_rate_metric =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "link_rate_Mbps", "Records per-link rate, only contains one value",
            "Link source", "Link destination");

static auto* bin_size_metric =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<uint64_t>(
        "bin_size_ms", "Records how long (in milliseconds) is each bin");

static auto* route_add_metric =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<uint64_t>(
        "route_add",
        "Records how many routes need to be added per optimization");

static auto* route_update_metric =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<uint64_t>(
        "route_update",
        "Records how many routes need to be updated per optimization");

static auto* route_remove_metric =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<uint64_t>(
        "route_remove",
        "Records how many routes need to be removed per optimization");

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
                 std::chrono::milliseconds total_duration,
                 size_t periods_in_history, RoutingSystem* routing_system)
    : periods_in_history_(periods_in_history),
      history_bin_size_(history_bin_size),
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
  bins_in_history_ = period_duration_bins_ * periods_in_history;
  CHECK(period_duration_bins_ > 0);
  period_count_ = min_bin_count / period_duration_bins_;

  period_count_ = std::min(
      period_count_,
      static_cast<size_t>(total_duration.count() / period_duration.count()));

  bin_size_metric->GetHandle()->AddValue(bin_size.count());
  for (nc::net::GraphLinkIndex link : graph_->AllLinks()) {
    const nc::net::GraphLink* link_ptr = graph_->GetLink(link);
    link_rate_metric->GetHandle(link_ptr->src_id(), link_ptr->dst_id())
        ->AddValue(link_ptr->bandwidth().Mbps());
  }
}

// Generates the input to the system.
std::map<AggregateId, AggregateHistory> NetMock::GenerateInput(
    const std::map<AggregateId, BinSequence>& period_sequences,
    PcapDataBinCache* cache) const {
  using namespace std::chrono;
  auto t1 = high_resolution_clock::now();
  std::map<AggregateId, AggregateHistory> input;
  for (const auto& aggregate_and_bins : period_sequences) {
    const AggregateId& aggregate = aggregate_and_bins.first;
    const BinSequence& bins = aggregate_and_bins.second;

    input.emplace(std::piecewise_construct, std::forward_as_tuple(aggregate),
                  std::forward_as_tuple(
                      bins.GenerateHistory(history_bin_size_, 1000, cache)));
  }
  auto t2 = high_resolution_clock::now();
  LOG(INFO) << "Gen input in " << duration_cast<milliseconds>(t2 - t1).count();

  return input;
}

nc::net::GraphLinkMap<std::vector<std::pair<double, double>>>
NetMock::CheckOutput(const std::map<AggregateId, BinSequence>& period_sequences,
                     const RoutingConfiguration& configuration,
                     PcapDataBinCache* cache) const {
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

  nc::net::GraphLinkMap<std::vector<std::pair<double, double>>> out;
  for (const auto& link_and_bins : link_to_bins) {
    nc::net::GraphLinkIndex link = link_and_bins.first;
    nc::net::Bandwidth rate = graph_->GetLink(link)->bandwidth();
    out[link] = (*link_and_bins.second)->Residuals(rate, cache);
  }

  return out;
}

std::map<AggregateId, BinSequence> NetMock::GetNthPeriod(
    size_t n, size_t bin_count) const {
  size_t period_start_bin = n * period_duration_bins_;
  size_t period_end_bin = period_start_bin + bin_count;

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
  std::map<AggregateId, BinSequence> zero_period =
      GetNthPeriod(0, bins_in_history_);
  std::map<AggregateId, AggregateHistory> input =
      GenerateInput(zero_period, cache);
  return routing_system_->Update(input).routing;
}

static size_t CheckSameSize(const nc::net::GraphLinkMap<
                            std::vector<std::pair<double, double>>>& values) {
  size_t i = 0;
  for (const auto& link_and_values : values) {
    const auto& v = *link_and_values.second;
    CHECK(v.size() > 0);
    if (i == 0) {
      i = v.size();
    } else {
      CHECK(i == v.size());
    }
  }

  return i;
}

std::map<AggregateId, nc::net::Bandwidth> NetMock::GetMeansForNthPeriod(
    size_t n, PcapDataBinCache* cache) const {
  std::map<AggregateId, BinSequence> period_sequences =
      GetNthPeriod(n, bins_in_history_);
  std::map<AggregateId, nc::net::Bandwidth> out;
  for (const auto& id_and_sequence : period_sequences) {
    const AggregateId& id = id_and_sequence.first;
    nc::net::Bandwidth mean_level = id_and_sequence.second.MeanRate(cache);
    if (mean_level == nc::net::Bandwidth::Zero()) {
      mean_level = nc::net::Bandwidth::FromBitsPerSecond(10);
    }

    out[id] = mean_level;
  }

  return out;
}

static void RecordDelta(uint64_t at, const RoutingConfiguration& from,
                        const RoutingConfiguration& to) {
  RoutingConfigurationDelta delta = from.GetDifference(to);

  uint64_t total_add = 0;
  uint64_t total_update = 0;
  uint64_t total_remove = 0;
  for (const auto& aggregate_and_delta : delta.aggregates) {
    const AggregateDelta& delta = aggregate_and_delta.second;
    total_add += delta.routes_added;
    total_update += delta.routes_updated;
    total_remove += delta.routes_removed;
  }

  route_add_metric->GetHandle()->AddValueWithTimestamp(at, total_add);
  route_update_metric->GetHandle()->AddValueWithTimestamp(at, total_update);
  route_remove_metric->GetHandle()->AddValueWithTimestamp(at, total_remove);
}

void NetMock::Run(PcapDataBinCache* cache) {
  std::unique_ptr<RoutingConfiguration> output = InitialOutput(cache);
  size_t timestamp = 0;
  for (size_t i = 0; i < period_count_; ++i) {
    LOG(ERROR) << "Period " << i;

    std::map<AggregateId, BinSequence> period_sequences =
        GetNthPeriod(i, period_duration_bins_);
    nc::net::GraphLinkMap<std::vector<std::pair<double, double>>>
        per_link_residuals = CheckOutput(period_sequences, *output, cache);
    size_t num_residuals = CheckSameSize(per_link_residuals);

    for (const auto& link_and_residuals : per_link_residuals) {
      nc::net::GraphLinkIndex link = link_and_residuals.first;
      const std::vector<std::pair<double, double>>& bins_and_residuals =
          *link_and_residuals.second;

      const nc::net::GraphLink* link_ptr = graph_->GetLink(link);
      auto* link_utilization_handle = link_utilization_metric->GetHandle(
          link_ptr->src_id(), link_ptr->dst_id());
      auto* link_queue_size_handle =
          queue_size_metric->GetHandle(link_ptr->src_id(), link_ptr->dst_id());

      size_t t = timestamp;
      for (size_t i = 0; i < bins_and_residuals.size(); ++i) {
        link_utilization_handle->AddValueWithTimestamp(
            t, bins_and_residuals[i].first);
        link_queue_size_handle->AddValueWithTimestamp(
            t++, bins_and_residuals[i].second);
      }
    }

    timestamp += num_residuals;

    // At the end of the period, will run the optimization or will check to see
    // if we need to trigger an optimization.
    size_t next_period = i + 1;
    if (next_period < periods_in_history_) {
      // Skip history for the first period.
      continue;
    }

    std::map<AggregateId, BinSequence> sequences_for_last_history =
        GetNthPeriod(next_period - periods_in_history_, bins_in_history_);
    std::map<AggregateId, AggregateHistory> input =
        GenerateInput(sequences_for_last_history, cache);

    bool need_update = false;
    if (next_period % periods_in_history_ == 0) {
      LOG(INFO) << "Will force update";
      need_update = true;
    } else if (FLAGS_triggered_optimization) {
      std::set<AggregateId> aggregates_no_fit;
      std::tie(aggregates_no_fit, std::ignore) =
          routing_system_->CheckWithProbModel(*output, input);
      if (!aggregates_no_fit.empty()) {
        LOG(INFO) << "Will trigger update";
        need_update = true;
      }
    }

    if (!need_update) {
      continue;
    }

    RoutingSystemUpdateResult result;
    if (FLAGS_mean_hints) {
      std::map<AggregateId, nc::net::Bandwidth> next_period_means =
          GetMeansForNthPeriod(next_period, cache);
      result = routing_system_->Update(input, next_period_means);
    } else {
      result = routing_system_->Update(input);
    }

    RecordDelta(timestamp, *output, *result.routing);
    output = std::move(result.routing);
  }
}

}  // namespace ctr
