#ifndef CTR_CTR_H
#define CTR_CTR_H

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include "ncode_net/src/net_common.h"
#include "../common.h"
#include "opt.h"

namespace ctr {

// Paths to be used by the optimizer. For each aggregate a list of paths in
// increasing order of delay.
using CTRPathMap = std::map<AggregateId, std::vector<const nc::net::Walk*>>;

struct RunOutput {
  RunOutput()
      : obj_value(std::numeric_limits<double>::max()),
        max_oversubscription(0) {}

  std::map<AggregateId, std::vector<RouteAndFraction>> aggregate_outputs;
  double obj_value;
  double max_oversubscription;
};

// A single pass of the optimizer.
class CTROptimizerPass {
 public:
  CTROptimizerPass(const TrafficMatrix* input, const CTRPathMap* paths,
                   const nc::net::GraphStorage* graph,
                   double link_capacity_multiplier);

  const nc::net::GraphLinkSet& links_with_no_capacity() const {
    return links_with_no_capacity_;
  }

  RunOutput& run_output() { return run_output_; }

 private:
  using FrozenCapacityMap = nc::net::GraphLinkMap<double>;

  // Freezes single path aggregates and adds them to the solution.
  void FreezeSinglePathAggregates();

  // Performs a single pass by repeatedly calling
  // OptimizerMinLinkOversubscription.
  void Optimize();

  // Runs the optimization with a given set of paths. Frozen aggregates will be
  // excluded from the optimization. Frozen per-link capacity will be excluded
  // from a link's capacity. Aggregates that have all of their paths go over
  // links at max oversubscription will be frozen after the optimization. The
  // freeze_all argument will cause all aggregates to be frozen after the
  // optimization, regardless of where their paths go.
  double OptimizeMinLinkOversubscription();

  // The input. Not owned by this class.
  const TrafficMatrix* input_;

  // A map of paths that the solver can use. Not owned by this class.
  const CTRPathMap* paths_;

  // The graph.
  const nc::net::GraphStorage* graph_;

  // Aggregates which are frozen are not going to be optimized by the solver.
  std::set<AggregateId> frozen_aggregates_;

  // Per-link capacity that the solver should ignore.
  FrozenCapacityMap frozen_capacity_;

  // Value of max oversubscription after the most recent optimization run.
  double latest_run_max_oversubscription_;

  // Links that have no capacity left on them.
  nc::net::GraphLinkSet links_with_no_capacity_;

  // Stores the output.
  RunOutput run_output_;

  // Oversubscription after all single-path aggregates have been frozen.
  double initial_oversubscription_;

  // Objective function value after all single-path aggregates have been frozen.
  double initial_obj_;

  // All links capacities' will be multiplied by this number.
  double link_capacity_multiplier_;
};

class CTROptimizer : public Optimizer {
 public:
  CTROptimizer(std::unique_ptr<PathProvider> path_provider,
               double link_capacity_multiplier = 1.0)
      : Optimizer(std::move(path_provider), link_capacity_multiplier) {}

  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override;

 private:
  // Single run of the optimization for given input.
  void OptimizePrivate(const TrafficMatrix& input,
                       const std::vector<AggregateId>& aggregates_ordered,
                       RoutingConfiguration* out);

  // Prioritizes the aggregates from an input. Aggregates that are more
  // important are at the start of the returned list.
  std::vector<AggregateId> PrioritizeAggregates(
      const TrafficMatrix& input) const;

  // For each aggregate will add all paths up to and including the lowest delay
  // path that has some free capacity. Returns whether it added any paths.
  bool AddFreePaths(const nc::net::GraphLinkSet& links_with_no_capacity,
                    const std::vector<AggregateId>& aggregates_ordered,
                    CTRPathMap* out);

  // If the total number of paths that will go to the optimizer exceeds this
  // limit only the next best free path will be added, as opposed to all k
  // shortest paths until the free path.
  size_t soft_path_limit_ = 100000;

  // If the total number of paths that will go to the optimizer exceeds this
  // limit no more paths will be added (each aggregate after that will only have
  // one path -- its shortest).
  size_t hard_path_limit_ = 200000;

  // Each aggregate will only have up to this many K shortest paths added to it.
  // After this limit is reached, paths will be added in a way that does not
  // ensure optimality of the solution.
  size_t per_aggregate_path_limit_ = 1000;
};

}  // namespace ctr

#endif
