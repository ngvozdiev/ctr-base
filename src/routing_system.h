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

class RoutingSystem {
 public:
  RoutingSystem(const ProbModelConfig& prob_model_config, Optimizer* optimizer,
                PerAggregateMeanEstimatorFactory* mean_estimator_factory,
                bool store_to_metrics = true)
      : prob_model_config_(prob_model_config),
        optimizer_(optimizer),
        estimator_(mean_estimator_factory),
        graph_(optimizer_->graph()),
        store_to_metrics_(store_to_metrics) {}

  virtual ~RoutingSystem() {}

  std::unique_ptr<RoutingConfiguration> Update(
      const std::map<AggregateId, AggregateHistory>& history);

  const nc::net::GraphStorage* graph() const { return graph_; }

 private:
  // Fraction by which to scale up aggregates that do not fit.
  static constexpr double kScaleFraction = 0.1;

  // Checks if all aggregates fit. Returns the set of aggregates that do not
  // fit.
  std::set<AggregateId> CheckWithProbModel(
      const RoutingConfiguration& routing,
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

  ProbModelConfig prob_model_config_;

  // Computes a new optimal solution.
  Optimizer* optimizer_;

  // Predict each aggregate's level.
  MeanEstimator estimator_;

  const nc::net::GraphStorage* graph_;

  // If true will record to metrics.
  bool store_to_metrics_;
};

}  // namespace ctr

#endif
