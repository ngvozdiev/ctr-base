#ifndef CTR_MEAN_EST_H
#define CTR_MEAN_EST_H

#include <stddef.h>
#include <map>
#include <memory>
#include <vector>

#include "../common.h"

namespace ctr {

class PerAggregateMeanEstimator {
 public:
  virtual ~PerAggregateMeanEstimator() {}
  virtual double MeanForNextTimestep(
      const std::vector<double>& prev_values,
      const std::vector<double>& prev_estimates) = 0;
};

class PerAggregateMeanEstimatorFactory {
 public:
  virtual ~PerAggregateMeanEstimatorFactory() {}
  virtual std::unique_ptr<PerAggregateMeanEstimator> NewEstimator(
      const AggregateId& aggregate) = 0;
};

// Comes up with a mean value for the coming minute, based on means from
// previous minutes.
class MeanEstimator {
 public:
  // Scales aggregates so that they have a new mean by shifting the bins up or
  // down.
  static std::map<AggregateId, AggregateHistory> ScaleMeans(
      const std::map<AggregateId, AggregateHistory>& current,
      const std::map<AggregateId, nc::net::Bandwidth>& new_means);

  MeanEstimator(PerAggregateMeanEstimatorFactory* estimator_factory)
      : estimator_factory_(estimator_factory) {}

  std::map<AggregateId, AggregateHistory> EstimateNext(
      const std::map<AggregateId, AggregateHistory>& current);

 private:
  struct AggregateState {
    std::vector<double> previous_values;
    std::vector<double> previous_estimates;
    std::unique_ptr<PerAggregateMeanEstimator> mean_estimator;
  };

  // State associated with each aggregate.
  std::map<AggregateId, AggregateState> aggregates_;

  // Produces estimators on demand.
  PerAggregateMeanEstimatorFactory* estimator_factory_;
};

class MeanScaleEstimatorConfig {
 public:
  MeanScaleEstimatorConfig(double fixed_headroom, double min_decay_multiplier,
                           double max_decay_multiplier, size_t window);

  MeanScaleEstimatorConfig() : MeanScaleEstimatorConfig(0, 0, 0, 0) {}

  double fixed_headroom() const { return fixed_headroom_; }

  double max_decay_multiplier() const { return max_decay_multiplier_; }

  size_t window() const { return window_; }

  double min_decay_multiplier() const { return min_decay_multiplier_; }

 private:
  double fixed_headroom_;
  size_t window_;
  double min_decay_multiplier_;
  double max_decay_multiplier_;
};

// Scales the vale of the previous value timestep a fixed amount.
class MeanScaleEstimator : public PerAggregateMeanEstimator {
 public:
  MeanScaleEstimator(const MeanScaleEstimatorConfig& config)
      : config_(config), level_(0) {}

  double MeanForNextTimestep(
      const std::vector<double>& prev_values,
      const std::vector<double>& prev_estimates) override;

  const std::vector<double>& decay_multipliers() const {
    return decay_multipliers_;
  }

 private:
  double RatioInWindow(const std::vector<double>& prev_values,
                       const std::vector<double>& prev_estimates);

  const MeanScaleEstimatorConfig config_;
  double level_;
  std::vector<double> decay_multipliers_;
};

class MeanScaleEstimatorFactory : public PerAggregateMeanEstimatorFactory {
 public:
  MeanScaleEstimatorFactory(const MeanScaleEstimatorConfig& config)
      : config_(config) {}

  virtual std::unique_ptr<PerAggregateMeanEstimator> NewEstimator(
      const AggregateId& aggregate) override {
    nc::Unused(aggregate);
    return nc::make_unique<MeanScaleEstimator>(config_);
  }

 private:
  const MeanScaleEstimatorConfig config_;
};

}  // namespace ctr

#endif
