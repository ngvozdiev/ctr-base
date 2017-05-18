#ifndef CTR_COMMON_H
#define CTR_COMMON_H

#include <string>
#include <tuple>

#include "ncode_common/src/common.h"
#include "ncode_net/src/net_common.h"

namespace ctr {

// Each aggregate is identified by a combination of src, dst.
class AggregateId {
 public:
  AggregateId(nc::net::GraphNodeIndex src, nc::net::GraphNodeIndex dst)
      : src_(src), dst_(dst) {}

  nc::net::GraphNodeIndex src() const { return src_; }

  nc::net::GraphNodeIndex dst() const { return dst_; }

  std::string ToString(const nc::net::GraphStorage& graph) const;

  friend bool operator<(const AggregateId& a, const AggregateId& b);
  friend bool operator==(const AggregateId& a, const AggregateId& b);
  friend bool operator!=(const AggregateId& a, const AggregateId& b);

 private:
  nc::net::GraphNodeIndex src_;
  nc::net::GraphNodeIndex dst_;
};

// A bandwidth demand and arbitrary priority.
using DemandAndPriority = std::pair<nc::net::Bandwidth, double>;

// A path and a fraction of demand that should go over it.
using RouteAndFraction = std::pair<const nc::net::Walk*, double>;

// For each aggregate, a bandwidth demand and a priority.
class TrafficMatrix {
 public:
  // Loads a TraffixMatrix from a string in the format used by
  // https://bitbucket.org/StevenGay/repetita/src.
  static std::unique_ptr<TrafficMatrix> LoadRepetitaOrDie(
      const std::string& matrix_string,
      const std::vector<std::string>& node_names,
      const nc::net::GraphStorage* graph);

  // Constructs an empty traffic matrix.
  TrafficMatrix(const nc::net::GraphStorage* graph) : graph_(graph) {}

  const std::map<AggregateId, DemandAndPriority>& demands() const {
    return demands_;
  }

  void AddDemand(const AggregateId& aggregate_id,
                 const DemandAndPriority& demand_and_priority);

  // A DemandMatrix is similar to TrafficMatrix, but has no concept of priority.
  std::unique_ptr<nc::lp::DemandMatrix> ToDemandMatrix() const;

  // Returns a new traffic matrix with the same aggregates, but each aggregate's
  // demand is +- a fraction of this one.
  std::unique_ptr<TrafficMatrix> Randomize(double fraction,
                                           std::mt19937* rnd) const;

 private:
  // The graph.
  const nc::net::GraphStorage* graph_;

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
      const std::vector<RouteAndFraction>& routes_and_fractions);

  const std::vector<RouteAndFraction>& FindRoutesOrDie(
      const AggregateId& aggregate_id) const;

  const std::map<AggregateId, std::vector<RouteAndFraction>>& routes() const {
    return configuration_;
  }

  std::string ToString(const nc::net::GraphStorage& graph) const;

 private:
  std::map<AggregateId, std::vector<RouteAndFraction>> configuration_;

  DISALLOW_COPY_AND_ASSIGN(RoutingConfiguration);
};

}  // namespace ctr
#endif
