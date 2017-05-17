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
  TrafficMatrix() {}

  const std::map<AggregateId, DemandAndPriority>& demands() const {
    return demands_;
  }

  void AddDemand(const AggregateId& aggregate_id,
                 const DemandAndPriority& demand_and_priority) {
    CHECK(!nc::ContainsKey(demands_, aggregate_id));
    demands_[aggregate_id] = demand_and_priority;
  }

 private:
  // For each aggregate its demand and its priority.
  std::map<AggregateId, DemandAndPriority> demands_;

  DISALLOW_COPY_AND_ASSIGN(TrafficMatrix);
};

// For each aggregate a set of paths and a fraction of demand to route.
class RoutingConfiguration {
 public:
  RoutingConfiguration() {}

  void AddRouteAndFraction(
      const AggregateId& aggregate_id,
      const std::vector<RouteAndFraction>& routes_and_fractions) {
    CHECK(!nc::ContainsKey(configuration_, aggregate_id));

    // Just to be on the safe side will check that the sum of all fractions is 1
    double total = 0.0;
    for (const auto& route_and_fraction : routes_and_fractions) {
      total += route_and_fraction.second;
    }
    CHECK(total <= 1.001 && total >= 0.999);
    configuration_[aggregate_id] = routes_and_fractions;
  }

  const std::vector<RouteAndFraction>& FindRoutesOrDie(
      const AggregateId& aggregate_id) const {
    return nc::FindOrDieNoPrint(configuration_, aggregate_id);
  }

  const std::map<AggregateId, std::vector<RouteAndFraction>>& routes() const {
    return configuration_;
  }

 private:
  std::map<AggregateId, std::vector<RouteAndFraction>> configuration_;

  DISALLOW_COPY_AND_ASSIGN(RoutingConfiguration);
};

}  // namespace ctr
#endif
