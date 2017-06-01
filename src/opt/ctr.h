#ifndef CTR_CTR_H
#define CTR_CTR_H

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include "ncode_common/src/net/net_common.h"
#include "../common.h"
#include "opt.h"
#include "oversubscription_model.h"

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
                   const RoutingConfiguration* base_solution);

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

  // Returns a limit on the fraction of an aggregate's traffic that can go on a
  // path. Will use base_solution_ to do so if available. If base_solution_ is
  // null will always return 1.
  double PathLimitFraction(const AggregateId& aggregate,
                           const nc::net::Walk* path) const;

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

  // Each aggregate's path fractions will be set to never exceed the ones in
  // this solution, if present. The exception is the aggregate's longest path in
  // base_solution_. Paths that are not in base_solution_ are not limited.
  const RoutingConfiguration* base_solution_;
};

// A heuristic that tries to generate a new solution from a previous one with as
// little change as possible.
class CTRQuickOptimizer {
 public:
  CTRQuickOptimizer(PathProvider* path_provider,
                    const RoutingConfiguration* previous)
      : path_provider_(path_provider),
        graph_(previous->graph()),
        previous_(previous),
        model_(*previous_) {}

  std::unique_ptr<RoutingConfiguration> Optimize(const TrafficMatrix& tm);

 private:
  using PathAndLoad = std::pair<const nc::net::Walk*, nc::net::Bandwidth>;

  nc::net::Bandwidth MinFreeCapacity(
      const nc::net::Links& links,
      const nc::net::GraphLinkMap<nc::net::Bandwidth>& extra_capacities,
      const nc::net::GraphLinkMap<nc::net::Bandwidth>& slack_capacities) const;

  // Tries to fit 'bw' into a series of paths. Will update both 'paths' with the
  // extra bandwidth per path and 'bw' to reflect how much of it managed to fit.
  void FindRoom(
      const nc::net::GraphLinkMap<nc::net::Bandwidth>& slack_capacities,
      nc::net::GraphLinkMap<nc::net::Bandwidth>* extra_capacities,
      std::vector<PathAndLoad>* paths, nc::net::Bandwidth* bw) const;

  // Returns the set of links that have no capacity left.
  nc::net::GraphLinkSet LinksWithNoCapacity(
      const nc::net::GraphLinkMap<nc::net::Bandwidth>& extra_capacities,
      const nc::net::GraphLinkMap<nc::net::Bandwidth>& slack_capacities) const;

  // Returns the free capacity on a link.
  nc::net::Bandwidth FreeCapacityOnLink(
      nc::net::GraphLinkIndex link,
      const nc::net::GraphLinkMap<nc::net::Bandwidth>& extra_capacities,
      const nc::net::GraphLinkMap<nc::net::Bandwidth>& slack_capacities) const;

  // Provides paths.
  PathProvider* path_provider_;

  // The graph.
  const nc::net::GraphStorage* graph_;

  const RoutingConfiguration* previous_;

  // Helps to figure out how much free capacity there is along paths.
  OverSubModel model_;
};

class CTROptimizer : public Optimizer {
 public:
  explicit CTROptimizer(PathProvider* path_provider)
      : Optimizer(path_provider) {}

  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override;

  std::unique_ptr<RoutingConfiguration> OptimizeWithPrevious(
      const TrafficMatrix& tm, const RoutingConfiguration& previous) override;

 private:
  // Single run of the optimization for given input.
  double OptimizePrivate(const TrafficMatrix& input,
                         const std::vector<AggregateId>& aggregates_ordered,
                         const RoutingConfiguration* base_solution,
                         RoutingConfiguration* out);

  // Prioritizes the aggregates from an input. Aggregates that are more
  // important are at the start of the returned list.
  std::vector<AggregateId> PrioritizeAggregates(
      const TrafficMatrix& input) const;

  // For each aggregate will add all paths up to and including the lowest delay
  // path that has some free capacity. Returns whether it added any paths.
  bool AddFreePaths(const nc::net::GraphLinkSet& links_with_no_capacity,
                    const std::vector<AggregateId>& aggregates_ordered,
                    const RoutingConfiguration* base_solution,
                    std::map<AggregateId, size_t>* ksp_indices,
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
