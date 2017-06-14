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
  virtual ~RoutingSystem() {}

  std::unique_ptr<RoutingConfiguration> Update(
      const std::map<AggregateId, AggregateHistory>& history);

 private:
  // Checks if all aggregates fit. Returns the set of aggregates that do not
  // fit.
  std::set<AggregateId> CheckWithProbModel(
      const RoutingConfiguration& routing,
      const std::map<AggregateId, AggregateHistory>& histories);

  bool ScaleUpAggregates(
      const std::set<AggregateId>& aggregates,
      const std::map<AggregateId, AggregateHistory>& histories,
      std::map<AggregateId, DemandAndFlowCount>* out);

  // Generates an initial input where all aggregates' rates are set to their
  // mean rate.
  std::map<AggregateId, DemandAndFlowCount> GetInitialInput(
      const std::map<AggregateId, AggregateHistory>& histories);

  // Fraction by which to scale up aggregates that do not fit.
  double scale_fraction_;

  ProbModelConfig prob_model_config_;

  // Computes a new optimal solution.
  std::unique_ptr<Optimizer> optimizer_;

  // Predict each aggregate's level.
  std::unique_ptr<MeanEstimator> estimator_;

  const nc::net::GraphStorage* graph_;
};

}  // namespace ctr

#endif
