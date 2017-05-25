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
  Optimizer(std::unique_ptr<PathProvider> path_provider)
      : path_provider_(std::move(path_provider)) {
    CHECK(path_provider_);
    graph_ = path_provider_->graph();
  }

  virtual ~Optimizer() {}

  virtual std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) = 0;

 protected:
  // Provides paths to the optimizer.
  std::unique_ptr<PathProvider> path_provider_;

  // The graph.
  const nc::net::GraphStorage* graph_;
};

// Only ever assigns aggregates to their single shortest path.
class ShortestPathOptimizer : public Optimizer {
 public:
  ShortestPathOptimizer(std::unique_ptr<PathProvider> path_provider)
      : Optimizer(std::move(path_provider)) {}

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
  MinMaxOptimizer(std::unique_ptr<PathProvider> path_provider)
      : Optimizer(std::move(path_provider)) {}

  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override;
};

// Runs a heuristic similar to that of B4. Each aggregate's "fair share" will be
// set to its flow count from the TM.
class B4Optimizer : public Optimizer {
 public:
  B4Optimizer(std::unique_ptr<PathProvider> path_provider)
      : Optimizer(std::move(path_provider)) {}

  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override;
};

}  // namespace ctr

#endif
