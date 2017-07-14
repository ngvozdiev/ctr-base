#ifndef TLDR_H
#define TLDR_H

#include <stddef.h>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/packet.h"
#include "ncode_common/src/net/net_common.h"
#include "common.h"
#include "prob_model/dist_model.h"

namespace ctr {

struct TLDRConfig {
  TLDRConfig(const nc::ThresholdEnforcerPolicy& threshold_enforcer_policy,
             const ProbModelConfig& prob_model_config,
             nc::net::IPAddress ip_src, nc::net::IPAddress ip_switch_dst,
             nc::net::IPAddress ip_controller_dst,
             std::chrono::milliseconds round_len,
             std::chrono::milliseconds switch_poll_period,
             size_t flow_count_sample_n,
             bool disable_fast_optimization_requests)
      : ip_src(ip_src),
        ip_switch_dst(ip_switch_dst),
        ip_controller_dst(ip_controller_dst),
        round_len(round_len),
        threshold_enforcer_policy(threshold_enforcer_policy),
        switch_poll_period(switch_poll_period),
        flow_count_sample_n(flow_count_sample_n),
        disable_fast_optimization_requests(disable_fast_optimization_requests),
        prob_model_config(prob_model_config) {}

  // All messages will have this address as source.
  nc::net::IPAddress ip_src;

  // Messages for the switch will have this destination address.
  nc::net::IPAddress ip_switch_dst;

  // Messages for the controller will have this destination address.
  nc::net::IPAddress ip_controller_dst;

  // How often to send LDR estimates to the controller.
  std::chrono::milliseconds round_len;

  // Configures the threshold enforcer that limits churn to the switch.
  nc::ThresholdEnforcerPolicy threshold_enforcer_policy;

  // How often to poll the switch for per-rule stats.
  std::chrono::milliseconds switch_poll_period;

  // One in N packets will be sampled to estimate flow counts at the switch.
  size_t flow_count_sample_n;

  // If true will disable optimization requests within the period.
  bool disable_fast_optimization_requests;

  ProbModelConfig prob_model_config;
};

// Bins data to later be fed to the LDR estimator.
class Binner {
 public:
  virtual ~Binner() {}

  // Returns the last N bins.
  std::vector<uint64_t> GetLastNBins(size_t n) const;

  // Returns all bins since (and including) 'from'.
  std::vector<uint64_t> GetBinsSince(nc::EventQueueTime from) const;

  std::vector<uint64_t> GetAllBins() const {
    return GetLastNBins(std::numeric_limits<size_t>::max());
  }

  size_t bin_count() const { return bins_.size(); }

 protected:
  struct Bin {
    Bin(nc::EventQueueTime start)
        : start(start), end(nc::EventQueueTime::MaxTime()), bytes(0) {}

    Bin(nc::EventQueueTime start, nc::EventQueueTime end, size_t bytes)
        : start(start), end(end), bytes(bytes) {}

    nc::EventQueueTime start;
    nc::EventQueueTime end;
    size_t bytes;
  };

  Binner(nc::EventQueue* event_queue) : event_queue_(event_queue) {}

  // Records for each bin the time the bin started and the number of bytes.
  std::vector<Bin> bins_;

  nc::EventQueue* event_queue_;
};

// Estimator that can have variable bin sizes, bins are added explicitly, each
// bin's duration is computed from the end of the previous bin.
class VariableBinBinner : public Binner {
 public:
  VariableBinBinner(nc::EventQueue* event_queue, bool cumulative_bins)
      : Binner(event_queue), last_bin_size_(0), cumulative_(cumulative_bins) {}

  // Adds a new bin.
  void AddBin(size_t bin_bytes);

 private:
  size_t last_bin_size_;
  bool cumulative_;
};

class TLDR : public ::nc::htsim::PacketHandler, public ::nc::EventConsumer {
 public:
  TLDR(const std::string& id, const TLDRConfig& config,
       nc::htsim::PacketHandler* to_switch,
       nc::htsim::PacketHandler* to_controller,
       const nc::net::GraphStorage* graph_storage, nc::EventQueue* event_queue);

  void HandlePacket(::nc::htsim::PacketPtr pkt) override;

  void HandleEvent() override;

 private:
  struct TLDRAggregateState {
    // Bins incoming traffic.
    std::unique_ptr<VariableBinBinner> binner;

    // To make things more stable flows are only counted once per period, and
    // the value is cached, then the counter reset. Every time a flow count is
    // needed the cached value is used.
    uint64_t cached_flow_count;

    // The most recent update sent from the central controller. This contains
    // per-path limits.
    std::unique_ptr<AggregateUpdateState> most_recent_update;

    // Enforces changes that are sent to the switch.
    std::unique_ptr<nc::ThresholdEnforcer<nc::htsim::PacketTag>>
        threshold_enforcer;

    // Stores the currently installed actions for this aggregate. This is used
    // to make sure we never delete actions -- any new updates that we  send to
    // the switch must not delete actions, as this will lose counter values.
    std::set<std::pair<nc::net::DevicePortNumber, nc::htsim::PacketTag>>
        actions_installed;
  };

  static constexpr size_t kWeightBase = 1000;

  void SendRequestToController(bool quick);

  // Updates the per-path bps limit metrics.
  void UpdateRateLimitMetrics(
      const std::map<AggregateId, AggregateUpdateState>& aggregates);

  void HandleUpdate(const TLDRUpdate& update);

  void UpdateTagToAggregateId(
      const std::map<AggregateId, AggregateUpdateState>& aggregates);

  void RepackPaths(const AggregateUpdateState& aggregate_update_state);

  void RepackPaths(
      const std::map<AggregateId, AggregateUpdateState>& aggregates);

  void UpdateAggregateState(
      const std::map<AggregateId, AggregateUpdateState>& aggregates,
      const std::map<AggregateId, AggregateHistory>&
          competing_aggregate_histories);

  // Returns the fractional change (in range 0-1) of a path's fraction from the
  // currently assigned fraction.
  double FindChange(nc::htsim::PacketTag tag, double new_fraction);

  // Returns the history seen for a given aggregate.
  AggregateHistory GetHistoryForAggregate(
      const TLDRAggregateState& aggregate_state) const;

  void HandleStatsReplyNoFlowCounts(const nc::htsim::SSCPStatsReply& reply);

  void HandleStatsReplyFlowCounts(const nc::htsim::SSCPStatsReply& reply);

  // Checks against competing aggregates whether the most recent history, when
  // combined with the histories of competing aggregates is likely to cause
  // congestion.
  bool LikelyToGoOverCapacity(
      const AggregateId& aggregate_id,
      const AggregateHistory& most_recent_history,
      const AggregateUpdateState& most_recent_update_state);

  // A copy of the initial config.
  const TLDRConfig tldr_config_;

  // Channel to the switch that this object manages.
  nc::htsim::PacketHandler* to_switch_;

  // Channel to the central controller.
  nc::htsim::PacketHandler* to_controller_;

  // Relates packet tags (which indicate paths) to aggregate ids. Updated
  // every time a new TLDRUpdate is received.
  std::map<nc::htsim::PacketTag, AggregateId> tag_to_aggregate_id_;

  // Internal per-aggregate state, indexed by aggregate id.
  std::map<AggregateId, TLDRAggregateState> id_to_aggregate_state_;

  // The histories of competing aggregates are stored here.
  std::map<AggregateId, AggregateHistory>
      most_recent_competing_aggregate_histories_;

  // The graph.
  const nc::net::GraphStorage* graph_;

  // How often to query the device.
  nc::EventQueueTime device_poll_period_;

  // How many times the device should be polled each round.
  uint64_t device_polls_in_round_;

  // How many times the device has been polled. Needed in order to figure out
  // when to send a request that resets flow counts.
  uint64_t device_poll_count_;

  // Generates round ids.
  uint64_t round_id_gen_;

  // Time the last triggered update happened. Used to avoid swamping the
  // controller with tons of triggered updates.
  nc::EventQueueTime last_trigered_update_;
};

}  // namespace tldr
#endif
