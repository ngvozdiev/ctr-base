#ifndef CTR_OPT_H
#define CTR_OPT_H

#include <algorithm>
#include <memory>

#include "ncode_common/src/logging.h"
#include "../common.h"
#include "path_provider.h"

namespace ctr {

class Optimizer {
 public:
  explicit Optimizer(PathProvider* path_provider)
      : path_provider_(path_provider) {
    CHECK(path_provider_);
    graph_ = path_provider_->graph();
  }

  virtual ~Optimizer() {}

  virtual std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) = 0;

  // Same as optimize, but will also make use of a previous configuration. Not
  // supported by all optimizers.
  virtual std::unique_ptr<RoutingConfiguration> OptimizeWithPrevious(
      const TrafficMatrix& tm, const RoutingConfiguration& previous) {
    nc::Unused(previous);
    return Optimize(tm);
  }

 protected:
  // Provides paths to the optimizer.
  PathProvider* path_provider_;

  // The graph.
  const nc::net::GraphStorage* graph_;
};

// Only ever assigns aggregates to their single shortest path.
class ShortestPathOptimizer : public Optimizer {
 public:
  explicit ShortestPathOptimizer(PathProvider* path_provider)
      : Optimizer(path_provider) {}

  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override;
};

// Returns a solution that minimizes the maximum utilization of any link in the
// network. The solution will be optimal, and the paths that it uses will not
// come from the path provider (i.e., they will not conform to any policy).
// Should be considered a lower bound for how low link utilization can get with
// a given traffic matrix.
class MinMaxOptimizer : public Optimizer {
 public:
  explicit MinMaxOptimizer(PathProvider* path_provider)
      : Optimizer(path_provider) {}

  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override;
};

// Runs a heuristic similar to that of B4. Each aggregate's "fair share" will be
// set to its flow count from the TM.
class B4Optimizer : public Optimizer {
 public:
  explicit B4Optimizer(PathProvider* path_provider)
      : Optimizer(path_provider) {}

  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override;
};

}  // namespace ctr

#endif
