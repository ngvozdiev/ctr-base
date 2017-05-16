#ifndef CTR_COMMON_H
#define CTR_COMMON_H

#include "ncode_net/src/net_common.h"

namespace ctr {

using AggregateId = std::pair<nc::net::GraphNodeIndex, nc::net::GraphNodeIndex>;

// A bandwidth demand and arbitrary priority.
using DemandAndPriority = std::pair<nc::net::Bandwidth, double>;

// A path and a fraction of demand that should go over it.
using RouteAndFraction = std::pair<const nc::net::Walk*, double>;

// For each aggregate, a bandwidth demand and a priority.
class TrafficMatrix {
 public:
  const std::map<AggregateId, DemandAndPriority>& demands() const {
    return demands_;
  }

 private:
  // For each aggregate its demand and its priority.
  std::map<AggregateId, DemandAndPriority> demands_;
};

// For each aggregate a set of paths and a fraction of demand to route.
class RoutingConfiguration {
 public:
  void AddRouteAndFraction(const AggregateId& aggregate_id,
                           const RouteAndFraction& route_and_fraction) {
    CHECK(!nc::ContainsKey(configuration_, aggregate_id));
    configuration_[aggregate_id] = route_and_fraction;
  }

  RouteAndFraction FindRouteOrDie(const AggregateId& aggregate_id) const {
    return nc::FindOrDieNoPrint(configuration_, aggregate_id);
  }

 private:
  std::map<AggregateId, RouteAndFraction> configuration_;
};

}  // namespace ctr
