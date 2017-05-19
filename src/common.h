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
  // Constructs an empty traffic matrix.
  TrafficMatrix(const nc::net::GraphStorage* graph) : graph_(graph) {}

  // Constructs a traffic matrix from a demand matrix. The demand matrix has no
  // flow counts, so they need to be supplied explicitly per src/dst pair (the
  // second argument). If there is no entry in the map for a src/dst pair its
  // flow count is assumed to be 1.
  TrafficMatrix(const nc::lp::DemandMatrix& demand_matrix,
                const std::map<nc::lp::SrcAndDst, double>& flow_counts = {});

  const std::map<AggregateId, DemandAndFlowCount>& demands() const {
    return demands_;
  }

  void AddDemand(const AggregateId& aggregate_id,
                 const DemandAndFlowCount& demand_and_flow_count);

  // A DemandMatrix is similar to TrafficMatrix, but has no flow counts.
  std::unique_ptr<nc::lp::DemandMatrix> ToDemandMatrix() const;

  // Returns a new traffic matrix with the same aggregates, but each aggregate's
  // demand/flow count is +- a fraction of the demand/flow count in this one.
  std::unique_ptr<TrafficMatrix> Randomize(double demand_fraction,
                                           double flow_count_fraction,
                                           std::mt19937* rnd) const;

  const nc::net::GraphStorage* graph() const { return graph_; }

 private:
  // The graph.
  const nc::net::GraphStorage* graph_;

  // For each aggregate its demand and its flow count.
  std::map<AggregateId, DemandAndFlowCount> demands_;

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
