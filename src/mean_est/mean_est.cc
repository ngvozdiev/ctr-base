#include "mean_est.h"

#include <algorithm>
#include <utility>

#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/net_common.h"

namespace ctr {

std::map<AggregateId, AggregateHistory> MeanEstimator::EstimateNext(
    const std::map<AggregateId, AggregateHistory>& current) {
  std::map<AggregateId, AggregateHistory> out;
  for (const auto& aggregate_and_history : current) {
    const AggregateId& aggregate_id = aggregate_and_history.first;
    const AggregateHistory& history = aggregate_and_history.second;

    double rate_mbps = history.mean_rate().Mbps();
    AggregateState& state = aggregates_[aggregate_id];
    if (state.previous_estimates.empty()) {
      // Will bootstrap with the first estimate as current.
      state.previous_estimates.emplace_back(rate_mbps);
    }

    state.previous_values.emplace_back(rate_mbps);
    CHECK(state.previous_estimates.size() == state.previous_values.size());

    if (!state.mean_estimator) {
      state.mean_estimator = estimator_factory_->NewEstimator(aggregate_id);
    }

    double next_mbps = state.mean_estimator->MeanForNextTimestep(
        state.previous_values, state.previous_estimates);
    CHECK(next_mbps >= 0);
    state.previous_estimates.emplace_back(next_mbps);

    // Need to scale the AggregateHistory so that its mean is not at
    // rate_mbps, but at next_mbps.
    double delta_mbps = next_mbps - rate_mbps;

    if (delta_mbps > 0) {
      out.emplace(std::piecewise_construct, std::make_tuple(aggregate_id),
                  std::make_tuple(history.AddRate(
                      nc::net::Bandwidth::FromMBitsPerSecond(delta_mbps))));
    } else if (delta_mbps < 0) {
      out.emplace(std::piecewise_construct, std::make_tuple(aggregate_id),
                  std::make_tuple(history.SubtractRate(
                      nc::net::Bandwidth::FromMBitsPerSecond(-delta_mbps))));
    }
  }
  return out;
}

MeanScaleEstimatorConfig::MeanScaleEstimatorConfig(double fixed_headroom,
                                                   double min_decay_multiplier,
                                                   double max_decay_multiplier,
                                                   size_t window)
    : fixed_headroom_(fixed_headroom),
      window_(window),
      min_decay_multiplier_(min_decay_multiplier),
      max_decay_multiplier_(max_decay_multiplier) {}

double MeanScaleEstimator::MeanForNextTimestep(
    const std::vector<double>& prev_values,
    const std::vector<double>& prev_estimates) {
  double fixed_headroom = config_.fixed_headroom();
  CHECK(!prev_values.empty());

  double ratio = RatioInWindow(prev_values, prev_estimates);
  double decay_multiplier =
      config_.min_decay_multiplier() +
      (config_.max_decay_multiplier() - config_.min_decay_multiplier()) * ratio;
  decay_multipliers_.emplace_back(decay_multiplier);

  double scaled_back = prev_values.back() * fixed_headroom;
  if (scaled_back > level_) {
    level_ = scaled_back;
  } else {
    level_ *= decay_multiplier;
    level_ = std::max(level_, scaled_back);
  }

  return level_;
}

double MeanScaleEstimator::RatioInWindow(
    const std::vector<double>& prev_values,
    const std::vector<double>& prev_estimates) {
  CHECK(prev_values.size() == prev_estimates.size());
  size_t size = prev_values.size();

  double total_prev_values = 0;
  double total_prev_estimates = 0;

  size_t lookback = std::min(config_.window(), size);
  for (size_t i = size - lookback; i < size; ++i) {
    double prev_value = prev_values[i];
    double prev_estimate = prev_estimates[i];

    //    CHECK(prev_value <= prev_estimate) << "at step " << i << " " <<
    //    prev_value
    //                                       << " vs " << prev_estimate;
    total_prev_values += prev_value;
    total_prev_estimates += prev_estimate;
  }

  return std::max(0.0, total_prev_values / total_prev_estimates);
}

}  // namespace ctr
