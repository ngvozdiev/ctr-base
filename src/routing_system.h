#ifndef CTR_ROUTING_SYS_H
#define CTR_ROUTING_SYS_H

#include <cstdbool>
#include <map>
#include <memory>
#include <set>

#include "common.h"
#include "mean_est/mean_est.h"
#include "opt/opt.h"
#include "prob_model/dist_model.h"

namespace ctr {

struct RoutingSystemConfig {
  // Enables pinning aggregates' levels to either their mean or max value.
  // Enabling one of those disables scaling of aggregates based on a probability
  // model. Cannot both be true.
  bool pin_max = false;
  bool pin_mean = false;

  // If true will store information about the routing system to metrics.
  bool store_to_metrics = true;

  ProbModelConfig prob_model_config;
};

class RoutingSystem {
 public:
  RoutingSystem(const RoutingSystemConfig& config, Optimizer* optimizer,
                PerAggregateMeanEstimatorFactory* mean_estimator_factory)
      : config_(config),
        optimizer_(optimizer),
        estimator_(mean_estimator_factory),
        graph_(optimizer_->graph()) {}

  virtual ~RoutingSystem() {}

  std::pair<std::unique_ptr<RoutingConfiguration>,
            std::unique_ptr<CompetingAggregates>>
  Update(const std::map<AggregateId, AggregateHistory>& history);

  const nc::net::GraphStorage* graph() const { return graph_; }

 private:
  // Fraction by which to scale up aggregates that do not fit.
  static constexpr double kScaleFraction = 0.1;

  // Checks if all aggregates fit. Returns the set of aggregates that do not
  // fit.
  std::pair<std::set<AggregateId>, std::unique_ptr<CompetingAggregates>>
  CheckWithProbModel(const RoutingConfiguration& routing,
                     const std::map<AggregateId, AggregateHistory>& histories);

  bool ScaleUpAggregates(
      const std::set<AggregateId>& aggregates,
      const std::map<AggregateId, AggregateHistory>& histories,
      std::map<AggregateId, DemandAndFlowCount>* out,
      std::map<AggregateId, double>* scale_fractions);

  // Generates an initial input where all aggregates' rates are set to their
  // mean rate.
  std::map<AggregateId, DemandAndFlowCount> GetInitialInput(
      const std::map<AggregateId, AggregateHistory>& histories);

  RoutingSystemConfig config_;

  // Computes a new optimal solution.
  Optimizer* optimizer_;

  // Predict each aggregate's level.
  MeanEstimator estimator_;

  const nc::net::GraphStorage* graph_;
};

}  // namespace ctr

#endif
