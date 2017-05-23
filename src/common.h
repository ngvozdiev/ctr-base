#ifndef CTR_COMMON_H
#define CTR_COMMON_H

#include <stddef.h>
#include <random>
#include <string>

#include "ncode_common/src/lp/mc_flow.h"
#include "ncode_common/src/common.h"
#include "ncode_common/src/net/net_common.h"

namespace nc {
namespace lp {
class DemandMatrix;
} /* namespace lp */
} /* namespace nc */

namespace ctr {

// Each aggregate is identified by a combination of src, dst.
class AggregateId {
 public:
  AggregateId(nc::net::GraphNodeIndex src, nc::net::GraphNodeIndex dst)
      : src_(src), dst_(dst) {}

  nc::net::GraphNodeIndex src() const { return src_; }

  nc::net::GraphNodeIndex dst() const { return dst_; }

  std::string ToString(const nc::net::GraphStorage& graph) const;

  // Returns the delay of the shortest path through the graph that this
  // aggregate can be routed on.
  nc::net::Delay GetSPDelay(const nc::net::GraphStorage& graph) const;

  friend bool operator<(const AggregateId& a, const AggregateId& b);
  friend bool operator==(const AggregateId& a, const AggregateId& b);
  friend bool operator!=(const AggregateId& a, const AggregateId& b);

 private:
  nc::net::GraphNodeIndex src_;
  nc::net::GraphNodeIndex dst_;
};

// A bandwidth demand and flow count.
using DemandAndFlowCount = std::pair<nc::net::Bandwidth, size_t>;

// A path and a fraction of demand that should go over it.
using RouteAndFraction = std::pair<const nc::net::Walk*, double>;

// For each aggregate, a bandwidth demand and a flow count.
class TrafficMatrix {
 public:
  // Constructs an empty traffic matrix, or optionally let the caller
  // pre-populate it with demands and flow counts.
  explicit TrafficMatrix(
      const nc::net::GraphStorage* graph,
      const std::map<AggregateId, DemandAndFlowCount>& demand_and_counts = {})
      : graph_(graph) {
    demands_ = demand_and_counts;
  }

  // Constructs a traffic matrix from a demand matrix. The demand matrix has no
  // flow counts, so they need to be supplied explicitly per src/dst pair (the
  // second argument). If there is no entry in the map for a src/dst pair its
  // flow count is assumed to be 1.
  explicit TrafficMatrix(
      const nc::lp::DemandMatrix& demand_matrix,
      const std::map<nc::lp::SrcAndDst, size_t>& flow_counts = {});

  const std::map<AggregateId, DemandAndFlowCount>& demands() const {
    return demands_;
  }

  void AddDemand(const AggregateId& aggregate_id,
                 const DemandAndFlowCount& demand_and_flow_count);

  // A DemandMatrix is similar to TrafficMatrix, but has no flow counts.
  std::unique_ptr<nc::lp::DemandMatrix> ToDemandMatrix() const;

  // Returns a new traffic matrix with the same aggregates, but
  // 'aggregate_count' aggregates have their demand/flow count +- a fraction of
  // the demand/flow count in this one. If flow_count_fraction is set to
  // std::numeric_limits<double>::max() each aggregate's flow count in the new
  // matrix will be set in a way that preserves the demand/flow_count ratio in
  // this matrix.
  std::unique_ptr<TrafficMatrix> Randomize(double demand_fraction,
                                           double flow_count_fraction,
                                           size_t aggregate_count,
                                           std::mt19937* rnd) const;

  const nc::net::GraphStorage* graph() const { return graph_; }

  // Dumps the entire TM to string.
  std::string ToString() const;

  // Prints a summary of the TM.
  std::string SummaryToString() const;

 protected:
  // The graph.
  const nc::net::GraphStorage* graph_;

 private:
  // For each aggregate its demand and its flow count.
  std::map<AggregateId, DemandAndFlowCount> demands_;

  DISALLOW_COPY_AND_ASSIGN(TrafficMatrix);
};

// The difference between the same aggregate in two different outputs
// (RoutingConfigurations).
struct AggregateDelta {
  // The fraction of the aggregate's total traffic that changed paths.
  double fraction_delta = 0.0;

  // For before and after the path stretch is computed as the sum Ps * f over
  // all of the aggregate's paths where Ps is the path's stretch (absolute
  // difference from shortest path) and f is the fraction of the traffic that
  // goes on the path. The difference between the before and after quantities is
  // returned. If this value is positive the aggregate's flows will experience
  // less delay.
  double path_stretch_gain = 0.0;

  // New routes added (f == 0 -> f != 0).
  size_t routes_added = 0;

  // Old routes deleted (f != 0 -> f == 0).
  size_t routes_removed = 0;

  // Routes updated (f != 0 -> f != 0)
  size_t routes_updated = 0;
};

struct RoutingConfigurationDelta {
  // Fraction of total flows that changed path. Each aggregate's flow count is
  // assumed to be the max of the flow count before and after.
  double total_flow_fraction_delta;

  // Fraction of total volume that changed path. Each aggregate's volume is
  // assumed to be the max of the volume before and after.
  double total_volume_fraction_delta;

  // Per-aggregate deltas.
  std::map<AggregateId, AggregateDelta> aggregates;

  // Returns the number of routed added, removed and updated
  std::tuple<size_t, size_t, size_t> TotalRoutes() const;
};

// Extends a TM with for each aggregate a set of paths and a fraction of demand
// to route.
class RoutingConfiguration : public TrafficMatrix {
 public:
  explicit RoutingConfiguration(const TrafficMatrix& base_matrix)
      : TrafficMatrix(base_matrix.graph(), base_matrix.demands()) {}

  void AddRouteAndFraction(
      const AggregateId& aggregate_id,
      const std::vector<RouteAndFraction>& routes_and_fractions);

  const std::vector<RouteAndFraction>& FindRoutesOrDie(
      const AggregateId& aggregate_id) const;

  const std::map<AggregateId, std::vector<RouteAndFraction>>& routes() const {
    return configuration_;
  }

  std::string ToString() const;

  // Computes the difference between this routing configuration and another.
  // Both should have the same aggregates.
  RoutingConfigurationDelta GetDifference(
      const RoutingConfiguration& other) const;

 private:
  std::map<AggregateId, std::vector<RouteAndFraction>> configuration_;
};

}  // namespace ctr
#endif
