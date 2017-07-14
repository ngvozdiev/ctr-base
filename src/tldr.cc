#include "tldr.h"

#include <type_traits>

#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/strutil.h"
#include "metrics/metrics.h"

namespace ctr {

static auto* kPathFractionMetric =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<double, std::string>(
        "tldr_path_fraction", "Fraction sent to the switch", "path");

static auto* kMeanRateMeasuredMetric =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "tldr_measured_rate_mean",
            "Mean rate (in Mbps) measured by TLDR and sent to the controller",
            "aggregate source", "aggregate destination");

static auto* kMaxRateMeasuredMetric =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "tldr_measured_rate_max",
            "Max rate (in Mbps) measured by TLDR and sent to the controller",
            "aggregate source", "aggregate destination");

static auto* kFlowCountMeasuredMetric =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string, std::string>(
            "tldr_measured_flow_count",
            "Flow count measured by TLDR and sent to the controller",
            "aggregate source", "aggregate destination");

std::vector<uint64_t> Binner::GetLastNBins(size_t n) const {
  std::vector<uint64_t> out;
  size_t from = bins_.size() > n ? bins_.size() - n : 0;

  for (size_t i = from; i < bins_.size(); ++i) {
    const Bin& bin = bins_[i];
    out.emplace_back(bin.bytes);
  }

  return out;
}

std::vector<uint64_t> Binner::GetBinsSince(nc::EventQueueTime from) const {
  std::vector<uint64_t> out;

  for (const Bin& bin : bins_) {
    if (bin.start < from) {
      continue;
    }

    out.emplace_back(bin.bytes);
  }

  return out;
}

void VariableBinBinner::AddBin(size_t bin_size) {
  size_t to_add = bin_size;
  if (cumulative_) {
    CHECK(last_bin_size_ <= to_add) << last_bin_size_ << " vs " << to_add;
    to_add -= last_bin_size_;
  }

  auto now = event_queue_->CurrentTime();
  if (!bins_.empty()) {
    bins_.back().end = now;
  }

  bins_.emplace_back(now);
  bins_.back().bytes = to_add;
  last_bin_size_ = bin_size;
}

static bool MessageContainsFlowCounts(
    const nc::htsim::SSCPStatsReply& message) {
  for (const auto& key_and_stats : message.stats()) {
    const std::vector<nc::htsim::ActionStats>& action_stats =
        key_and_stats.second;

    for (const auto& stats : action_stats) {
      bool flow_count_meaningful =
          stats.flow_count != std::numeric_limits<uint64_t>::max();
      if (flow_count_meaningful) {
        return true;
      }
    }
  }

  return false;
}

bool TLDR::LikelyToGoOverCapacity(
    const AggregateId& aggregate_id,
    const AggregateHistory& most_recent_history,
    const AggregateUpdateState& most_recent_update_state) {
  if (most_recent_update_state.competing_aggregates.empty()) {
    return false;
  }

  ProbModel prob_model(tldr_config_.prob_model_config);
  prob_model.AddAggregate(aggregate_id, &most_recent_history);
  for (const auto& id_and_history :
       most_recent_competing_aggregate_histories_) {
    prob_model.AddAggregate(id_and_history.first, &id_and_history.second);
  }

  std::vector<ProbModelQuery> queries;
  for (size_t i = 0; i < most_recent_update_state.paths.size(); ++i) {
    const PathUpdateState& path_update_state =
        most_recent_update_state.paths[i];
    const AggregatesAndCapacity& competing_aggregates =
        most_recent_update_state.competing_aggregates[i];
    double fraction_on_path = path_update_state.route_and_fraction.second;

    queries.emplace_back();
    ProbModelQuery& query = queries.back();
    query.rate = competing_aggregates.capacity;
    query.type = ProbModelQuery::BOTH;
    query.aggregates = competing_aggregates.aggregates;
    query.aggregates.emplace_back(aggregate_id, fraction_on_path);
  }

  std::vector<ProbModelReply> replies = prob_model.Query(queries);
  for (const auto& reply : replies) {
    if (!reply.fit) {
      return true;
    }
  }

  return false;
}

void TLDR::HandleStatsReplyNoFlowCounts(
    const nc::htsim::SSCPStatsReply& stats_message) {
  std::map<AggregateId, uint64_t> id_to_bytes_matched;
  for (const auto& key_and_actions : stats_message.stats()) {
    for (const nc::htsim::ActionStats action_stats : key_and_actions.second) {
      // Need to figure out which aggregate the action belongs to by looking
      // at the tag.
      nc::htsim::PacketTag tag = action_stats.tag;
      AggregateId* id = nc::FindOrNull(tag_to_aggregate_id_, tag);
      if (id == nullptr) {
        continue;
      }

      id_to_bytes_matched[*id] += action_stats.total_bytes_matched;
    }
  }

  bool need_update = false;
  for (const auto& id_and_bytes_matched : id_to_bytes_matched) {
    AggregateId aggregate_id = id_and_bytes_matched.first;
    uint64_t bytes_matched = id_and_bytes_matched.second;

    TLDRAggregateState& aggregate_state =
        nc::FindOrDieNoPrint(id_to_aggregate_state_, aggregate_id);
    aggregate_state.binner->AddBin(bytes_matched);

    if (tldr_config_.disable_fast_optimization_requests) {
      continue;
    }

    // Now that we have updated the aggregate's rates we can check to see if
    // we need to send a message to the controller to force a new
    // optimization pass.
    AggregateHistory history = GetHistoryForAggregate(aggregate_state);
    if (history.bins().empty()) {
      continue;
    }

    CHECK(aggregate_state.most_recent_update);
    if (LikelyToGoOverCapacity(aggregate_id, history,
                               *aggregate_state.most_recent_update)) {
      need_update = true;
      LOG(INFO) << "Need update at " << id() << " aggregate id "
                << aggregate_id.ToString(*graph_);

      break;
    }
  }

  nc::EventQueueTime now = event_queue()->CurrentTime();
  CHECK(last_trigered_update_ < now);
  nc::EventQueueTime delta = now - last_trigered_update_;
  if (delta > event_queue()->ToTime(std::chrono::milliseconds(100)) &&
      need_update) {
    LOG(ERROR) << "Will trigger update at " << id();
    auto msg = nc::GetFreeList<TLDRTriggerReoptimize>().New(
        tldr_config_.ip_src, tldr_config_.ip_controller_dst,
        event_queue()->CurrentTime());
    to_controller_->HandlePacket(std::move(msg));
    last_trigered_update_ = event_queue()->CurrentTime();
  }
}

void TLDR::HandleStatsReplyFlowCounts(
    const nc::htsim::SSCPStatsReply& stats_message) {
  std::map<AggregateId, std::pair<uint64_t, uint64_t>>
      id_to_bytes_and_flow_counts;
  for (const auto& key_and_actions : stats_message.stats()) {
    for (const nc::htsim::ActionStats action_stats : key_and_actions.second) {
      nc::htsim::PacketTag tag = action_stats.tag;
      AggregateId* id = nc::FindOrNull(tag_to_aggregate_id_, tag);
      if (id == nullptr) {
        continue;
      }

      std::pair<uint64_t, uint64_t>& bytes_and_flow_counts =
          id_to_bytes_and_flow_counts[*id];

      bytes_and_flow_counts.first += action_stats.total_bytes_matched;
      bytes_and_flow_counts.second += action_stats.flow_count;
    }
  }

  for (const auto& id_and_bytes_and_flow_counts : id_to_bytes_and_flow_counts) {
    AggregateId aggregate_id = id_and_bytes_and_flow_counts.first;
    uint64_t bytes_matched, flow_count;
    std::tie(bytes_matched, flow_count) = id_and_bytes_and_flow_counts.second;

    TLDRAggregateState& aggregate_state =
        nc::FindOrDieNoPrint(id_to_aggregate_state_, aggregate_id);
    aggregate_state.binner->AddBin(bytes_matched);
    aggregate_state.cached_flow_count = flow_count;
  }

  SendRequestToController(false);
}

void TLDR::HandlePacket(::nc::htsim::PacketPtr pkt) {
  using nc::htsim::SSCPAddOrUpdate;
  using nc::htsim::SSCPAck;
  using nc::htsim::SSCPStatsReply;

  CHECK(pkt->size_bytes() == 0) << "Non-message packet at TLDR";

  uint8_t type = pkt->five_tuple().ip_proto().Raw();
  if (type == TLDRUpdate::kTLDRUpdateType) {
    TLDRUpdate* update_message = static_cast<TLDRUpdate*>(pkt.get());
    HandleUpdate(*update_message);
  } else if (type == SSCPAddOrUpdate::kSSCPAddOrUpdateType) {
    // SSPC messages / replies are directly forwarded to the switch /
    // controller.
    SSCPAddOrUpdate* original_update = static_cast<SSCPAddOrUpdate*>(pkt.get());

    auto update = nc::GetFreeList<SSCPAddOrUpdate>().New(
        tldr_config_.ip_src, tldr_config_.ip_switch_dst, pkt->time_sent(),
        original_update->TakeRule());
    update->set_tx_id(original_update->tx_id());
    to_switch_->HandlePacket(std::move(update));
  } else if (type == SSCPAck::kSSCPAckType) {
    SSCPAck* original_ack = static_cast<SSCPAck*>(pkt.get());
    auto ack = nc::GetFreeList<SSCPAck>().New(
        tldr_config_.ip_src, tldr_config_.ip_controller_dst, pkt->time_sent(),
        original_ack->tx_id());
    to_controller_->HandlePacket(std::move(ack));
  } else if (type == SSCPStatsReply::kSSCPStatsReplyType) {
    SSCPStatsReply* stats_message = static_cast<SSCPStatsReply*>(pkt.get());

    if (MessageContainsFlowCounts(*stats_message)) {
      HandleStatsReplyFlowCounts(*stats_message);
    } else {
      HandleStatsReplyNoFlowCounts(*stats_message);
    }

  } else if (type == TLDRForceRequest::kTLDRForceRequestType) {
    SendRequestToController(true);
  }
}

AggregateHistory TLDR::GetHistoryForAggregate(
    const TLDRAggregateState& aggregate_state) const {
  auto now = event_queue()->CurrentTime();
  nc::EventQueueTime time_for_round =
      device_poll_period_ * device_polls_in_round_;
  if (now < time_for_round) {
    return {{}, std::chrono::milliseconds::max(), 0};
  }
  std::vector<uint64_t> bins =
      aggregate_state.binner->GetBinsSince(now - time_for_round);
  if (bins.empty()) {
    return {{}, std::chrono::milliseconds::max(), 0};
  }

  uint64_t flow_count = aggregate_state.cached_flow_count;
  return {bins, tldr_config_.switch_poll_period, flow_count};
}

void TLDR::SendRequestToController(bool quick) {
  std::map<AggregateId, AggregateHistory> aggregate_to_history;
  for (auto& id_and_estimator : id_to_aggregate_state_) {
    AggregateId id = id_and_estimator.first;
    TLDRAggregateState& tldr_aggregate_state = id_and_estimator.second;
    AggregateHistory history = GetHistoryForAggregate(tldr_aggregate_state);
    if (history.bins().empty()) {
      continue;
    }

    std::string src = graph_->GetNode(id.src())->id();
    std::string dst = graph_->GetNode(id.dst())->id();

    aggregate_to_history.emplace(id, history);
    kFlowCountMeasuredMetric->GetHandle(src, dst)
        ->AddValue(history.flow_count());
    kMeanRateMeasuredMetric->GetHandle(src, dst)
        ->AddValue(history.mean_rate().Mbps());
    kMaxRateMeasuredMetric->GetHandle(src, dst)
        ->AddValue(history.max_rate().Mbps());
  }

  if (!aggregate_to_history.empty()) {
    auto request = nc::GetFreeList<TLDRRequest>().New(
        tldr_config_.ip_src, tldr_config_.ip_controller_dst,
        event_queue()->CurrentTime(), ++round_id_gen_, quick,
        aggregate_to_history);
    to_controller_->HandlePacket(std::move(request));
  }
}

void TLDR::HandleEvent() {
  bool include_flow_counts = ++device_poll_count_ % device_polls_in_round_ == 0;
  auto request = nc::GetFreeList<nc::htsim::SSCPStatsRequest>().New(
      tldr_config_.ip_src, tldr_config_.ip_switch_dst,
      event_queue()->CurrentTime(), include_flow_counts);
  to_switch_->HandlePacket(std::move(request));
  EnqueueIn(device_poll_period_);
}

void TLDR::HandleUpdate(const TLDRUpdate& update) {
  LOG(INFO) << "Rx " << update.ToString();
  UpdateAggregateState(update.aggregates(), update.additional_histories());
  UpdateTagToAggregateId(update.aggregates());
  RepackPaths(update.aggregates());
  UpdateRateLimitMetrics(update.aggregates());
}

void TLDR::UpdateTagToAggregateId(
    const std::map<AggregateId, AggregateUpdateState>& aggregates) {
  for (const auto& id_and_update_state : aggregates) {
    const AggregateId& id = id_and_update_state.first;
    const AggregateUpdateState& update_state = id_and_update_state.second;

    for (const PathUpdateState& path : update_state.paths) {
      AggregateId* current_id = nc::FindOrNull(tag_to_aggregate_id_, path.tag);
      if (current_id != nullptr) {
        CHECK(id == *current_id) << "Id/tag mismatch: " << id.ToString(*graph_)
                                 << " vs " << current_id->ToString(*graph_)
                                 << " for path tag " << path.tag;
      }

      nc::htsim::PacketTag tag(path.tag);
      tag_to_aggregate_id_.emplace(tag, id);
    }
  }
}

void TLDR::UpdateRateLimitMetrics(
    const std::map<AggregateId, AggregateUpdateState>& aggregates) {
  for (const auto& id_and_update_state : aggregates) {
    const AggregateId& id = id_and_update_state.first;
    const AggregateUpdateState& update_state = id_and_update_state.second;

    std::string src = graph_->GetNode(id.src())->id();
    std::string dst = graph_->GetNode(id.dst())->id();

    for (const PathUpdateState& path : update_state.paths) {
      const nc::net::Walk* graph_path = path.route_and_fraction.first;
      double fraction = path.route_and_fraction.second;

      std::string path_str = graph_path->ToStringNoPorts(*graph_);
      LOG(INFO) << path_str << " fraction " << fraction;
      kPathFractionMetric->GetHandle(path_str)->AddValue(fraction);
    }
  }
}

static void AddDummyActionToRule(nc::net::DevicePortNumber port,
                                 nc::htsim::PacketTag tag,
                                 nc::htsim::MatchRule* out) {
  for (const nc::htsim::MatchRuleAction* action : out->actions()) {
    if (action->tag() == tag && action->output_port() == port) {
      return;
    }
  }

  auto action = nc::make_unique<nc::htsim::MatchRuleAction>(port, tag, 0);
  out->AddAction(std::move(action));
}

void TLDR::RepackPaths(const AggregateUpdateState& aggregate) {
  TLDRAggregateState& current_state =
      nc::FindOrDieNoPrint(id_to_aggregate_state_, aggregate.aggregate_id);

  std::map<nc::htsim::PacketTag, double> tag_to_fraction;
  for (const PathUpdateState& path_state : aggregate.paths) {
    tag_to_fraction[path_state.tag] = path_state.route_and_fraction.second;
  }

  if (!current_state.threshold_enforcer->ChangeBulk(tag_to_fraction)) {
    // No significant change in any of the fractions, can ignore the update.
    LOG(INFO) << "Ignored change to aggregate "
              << aggregate.aggregate_id.ToString(*graph_);
    return;
  }

  auto match_rule = nc::make_unique<nc::htsim::MatchRule>(aggregate.key);
  CHECK(!aggregate.paths.empty());
  for (const PathUpdateState& path : aggregate.paths) {
    double fraction = nc::FindOrDie(tag_to_fraction, path.tag);

    size_t weight = kWeightBase * fraction;
    nc::htsim::PacketTag tag = path.tag;

    auto action = nc::make_unique<nc::htsim::MatchRuleAction>(
        path.output_port_num, tag, weight);
    current_state.actions_installed.emplace(path.output_port_num, tag);

    action->EnableFlowCounter(tldr_config_.flow_count_sample_n, event_queue());
    match_rule->AddAction(std::move(action));
  }

  // Need to populate the update with dummy actions for the ones that it is
  // missing. This will prevent counter values from being lost at the switch.
  for (const auto& output_port_and_tag : current_state.actions_installed) {
    AddDummyActionToRule(output_port_and_tag.first, output_port_and_tag.second,
                         match_rule.get());
  }

  auto route_update = nc::GetFreeList<nc::htsim::SSCPAddOrUpdate>().New(
      tldr_config_.ip_src, tldr_config_.ip_switch_dst,
      event_queue()->CurrentTime(), std::move(match_rule));
  to_switch_->HandlePacket(std::move(route_update));
}

void TLDR::RepackPaths(
    const std::map<AggregateId, AggregateUpdateState>& aggregates) {
  for (const auto& id_and_aggregate : aggregates) {
    RepackPaths(id_and_aggregate.second);
  }
}

void TLDR::UpdateAggregateState(
    const std::map<AggregateId, AggregateUpdateState>& aggregates,
    const std::map<AggregateId, AggregateHistory>&
        competing_aggregate_histories) {
  std::set<AggregateId> ids_in_update;
  for (const auto& id_and_update_state : aggregates) {
    const AggregateId& id = id_and_update_state.first;
    ids_in_update.insert(id);

    TLDRAggregateState* state_to_update;
    auto it = id_to_aggregate_state_.find(id);
    if (it != id_to_aggregate_state_.end()) {
      state_to_update = &(it->second);
    } else {
      state_to_update = &id_to_aggregate_state_[id];
      state_to_update->threshold_enforcer =
          nc::make_unique<nc::ThresholdEnforcer<nc::htsim::PacketTag>>(
              tldr_config_.threshold_enforcer_policy);
      state_to_update->binner =
          nc::make_unique<VariableBinBinner>(event_queue(), true);
      state_to_update->cached_flow_count = 1;
      LOG(INFO) << "Added state for aggregate " << id.ToString(*graph_);
    }

    state_to_update->most_recent_update =
        nc::make_unique<AggregateUpdateState>(id_and_update_state.second);
  }

  for (const auto& id_and_history : competing_aggregate_histories) {
    most_recent_competing_aggregate_histories_.emplace(id_and_history.first,
                                                       id_and_history.second);
  }

  for (auto it = id_to_aggregate_state_.begin();
       it != id_to_aggregate_state_.end(); ++it) {
    const AggregateId& aggregate_id = it->first;
    if (!nc::ContainsKey(ids_in_update, aggregate_id)) {
      id_to_aggregate_state_.erase(it);
      LOG(INFO) << "Removed state for aggregate "
                << aggregate_id.ToString(*graph_) << " at " << id();
      break;
    }
  }
}

TLDR::TLDR(const std::string& id, const TLDRConfig& config,
           nc::htsim::PacketHandler* to_switch,
           nc::htsim::PacketHandler* to_controller,
           const nc::net::GraphStorage* graph_storage,
           nc::EventQueue* event_queue)
    : EventConsumer(id, event_queue),
      tldr_config_(config),
      to_switch_(to_switch),
      to_controller_(to_controller),
      graph_(graph_storage),
      device_poll_period_(event_queue->ToTime(tldr_config_.switch_poll_period)),
      device_poll_count_(0),
      round_id_gen_(0),
      last_trigered_update_(nc::EventQueueTime::ZeroTime()) {
  size_t round_count = tldr_config_.round_len.count();
  size_t poll_period_count = tldr_config_.switch_poll_period.count();
  device_polls_in_round_ = round_count / poll_period_count;
  CHECK(round_count % poll_period_count == 0);
  EnqueueIn(device_poll_period_);
}
}  // namespace tldr
