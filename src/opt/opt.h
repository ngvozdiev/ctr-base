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

  const nc::net::GraphStorage* graph() { return path_provider_->graph(); }

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
  explicit B4Optimizer(PathProvider* path_provider,
                       bool flow_count_as_fair_share)
      : Optimizer(path_provider),
        flow_count_as_fair_share_(flow_count_as_fair_share) {}

  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override;

 private:
  bool flow_count_as_fair_share_;
};

}  // namespace ctr

#endif
