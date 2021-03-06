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
  // If true will store information about the routing system to metrics.
  bool store_to_metrics = true;

  // When performing variability tests each link's rate will be multiplied by
  // this number.
  double link_capacity_multiplier = 1.0;

  ProbModelConfig prob_model_config;
};

// The result of an update operation.
struct RoutingSystemUpdateResult {
  // The new routing configuration.
  std::unique_ptr<RoutingConfiguration> routing;

  // For each aggregate, the set of other aggregates it competes with and what
  // capacity it competes for.
  std::unique_ptr<CompetingAggregates> competing_aggregates;

  // The per-aggregate histories after prediction. Those will not be the same as
  // the histories before prediction (the input to the update operation).
  std::map<AggregateId, AggregateHistory> histories_after_prediction;
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

  RoutingSystemUpdateResult Update(
      const std::map<AggregateId, AggregateHistory>& history,
      const std::map<AggregateId, nc::net::Bandwidth>& mean_hints = {});

  const nc::net::GraphStorage* graph() const { return graph_; }

  // Checks if all aggregates fit. Returns the set of aggregates that do not
  // fit.
  std::pair<std::set<AggregateId>, std::unique_ptr<CompetingAggregates>>
  CheckWithProbModel(
      const RoutingConfiguration& routing,
      const std::map<AggregateId, AggregateHistory>& histories) const;

 private:
  // Fraction by which to scale up aggregates that do not fit.
  static constexpr double kScaleFraction = 0.1;

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
