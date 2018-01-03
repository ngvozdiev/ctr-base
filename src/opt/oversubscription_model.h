#ifndef FUBAR_OVERSUBSCRIPTION_MODEL_H
#define FUBAR_OVERSUBSCRIPTION_MODEL_H

#include <map>
#include <vector>

#include "ncode/net/net_common.h"
#include "../common.h"

namespace ctr {

class OverSubModelPathState;

// Helper per-link class used by OverSubModel.
class OverSubModelLinkState {
 public:
  OverSubModelLinkState(nc::net::Bandwidth capacity,
                        const nc::net::GraphLink* link);

  double SubscriptionRatio() const;

  void Bottleneck();

  void ReduceCapacity(nc::net::Bandwidth by) { remaining_capacity_ -= by; }

  void AddPathState(OverSubModelPathState* path_state);

  nc::net::Bandwidth remaining_capacity() const { return remaining_capacity_; }

  const std::vector<OverSubModelPathState*>& path_states_over_link() const {
    return path_states_over_link_;
  }

  const nc::net::GraphLink* link() const { return link_; }

  nc::net::Bandwidth original_capacity() const { return original_capacity_; }

 private:
  // The original capacity of this link state.
  const nc::net::Bandwidth original_capacity_;

  // The link that this state is for.
  const nc::net::GraphLink* link_;

  bool bottlenecked_;
  nc::net::Bandwidth remaining_capacity_;
  std::vector<OverSubModelPathState*> path_states_over_link_;
};

// Helper per-path class used by OverSubModel.
class OverSubModelPathState {
 public:
  OverSubModelPathState(const AggregateId& aggregate_id, double fraction,
                        nc::net::Bandwidth path_demand,
                        const nc::net::Walk* path);

  void Bottleneck(nc::net::Bandwidth at_rate);

  nc::net::Bandwidth bottleneck_rate() const;

  nc::net::Bandwidth path_demand() const { return path_demand_; }

  bool bottlenecked() const { return bottlenecked_; }

  const nc::net::Walk* path() const { return path_; }

  void AddLinkState(OverSubModelLinkState* link_state);

  const AggregateId& aggregate_id() const { return aggregate_id_; }

  double fraction() const { return fraction_; }

 private:
  double fraction_;

  // The id of the aggregate.
  AggregateId aggregate_id_;

  // The path this state is for.
  const nc::net::Walk* path_;

  bool bottlenecked_;
  nc::net::Bandwidth bottleneck_rate_;
  nc::net::Bandwidth path_demand_;
  std::vector<OverSubModelLinkState*> link_states_in_path_;
};

// When the topology cannot meet the demand of the TM some aggregates' demands
// will not be met. This class can be used to get a rough estimate for how
// loaded each link is and what its "oversubscription" value is.
class OverSubModel {
 public:
  OverSubModel(const RoutingConfiguration& routing);

  // The bandwidth a flow that is on a path is expected to get. If the network
  // is oversubscribed this will be less than the flow's demand (the total
  // demand of the aggregate divided by the number of flows in it).
  const std::map<const nc::net::Walk*, nc::net::Bandwidth>&
  per_flow_bandwidth_map() const {
    return per_flow_bandwidth_;
  }

  // Returns a map with each link's load. The load will never exceed 1.
  const nc::net::GraphLinkMap<double>& link_to_load() const {
    return link_to_load_;
  }

  // The set of aggregates that do not fit.
  const std::set<AggregateId>& aggregates_no_fit() const {
    return aggregates_no_fit_;
  }

 private:
  // Initializes state for a subset of the aggregates.
  void InitState(
      const RoutingConfiguration& routing,
      std::map<const nc::net::Walk*, OverSubModelPathState>* path_to_state,
      std::map<nc::net::GraphLinkIndex, OverSubModelLinkState>* link_to_state);

  // Will repeatedly find the link that is the most oversubscribed and
  // bottleneck it.
  void PopulateCapacities(
      std::map<const nc::net::Walk*, OverSubModelPathState>* path_to_state,
      std::map<nc::net::GraphLinkIndex, OverSubModelLinkState>* link_to_state);

  std::map<const nc::net::Walk*, nc::net::Bandwidth> per_flow_bandwidth_;

  nc::net::GraphLinkMap<double> link_to_load_;

  // The set of aggregates that have at least one path go over a link that does
  // not fit.
  std::set<AggregateId> aggregates_no_fit_;

  // The graph.
  const nc::net::GraphStorage* graph_;
};

}  // namespace fubarizer
#endif
