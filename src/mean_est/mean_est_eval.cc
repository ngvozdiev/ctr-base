#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <iterator>
#include <string>
#include <vector>

#include "ncode/common.h"
#include "ncode/file.h"
#include "ncode/logging.h"
#include "ncode/strutil.h"
#include "ncode/substitute.h"
#include "ncode/viz/grapher.h"
#include "mean_est.h"

DEFINE_string(time_series_file, "",
              "A comma-separated file with a single time series per line");
DEFINE_bool(
    cumulative, true,
    "Whether unmet demand will be retransmitted on the following timestep");
DEFINE_double(decay_fraction, 0.995, "How quickly the level decays");
DEFINE_double(fixed_scale, 1.1, "Fixed scale applied to all estimates");
DEFINE_uint64(window_size, 60, "How many timesteps in the past to examine");
DEFINE_bool(simple, true,
            "If true will run a simple estimator with fixed decay fraction");

namespace ctr {

class BaselineEstimator : public PerAggregateMeanEstimator {
 public:
  BaselineEstimator(const std::vector<double>& values) {
    CHECK(!values.empty());
    value_ = *std::max_element(values.begin(), values.end());
  }

  double MeanForNextTimestep(
      const std::vector<double>& prev_values,
      const std::vector<double>& prev_estimates) override {
    nc::Unused(prev_values);
    nc::Unused(prev_estimates);
    return value_;
  }

 private:
  double value_;
};
}

struct SimulationStats {
  double total_demand = 0;
  double total_unmet_demand = 0;
  double wasted_capacity = 0;
  size_t timesteps_with_unment_demand = 0;

  // Not really a statistic, but useful to have here. This is the estimator's
  // estimate for the future value.
  double next_value_estimate;

  std::string ToString() const {
    return nc::Substitute(
        "demand: $3, unmet: $0, wasted: $1, timesteps with unmet: $2",
        total_unmet_demand / total_demand, wasted_capacity / total_demand,
        timesteps_with_unment_demand, total_demand);
  }
};

// Simulates a process with per-timestep demand.
SimulationStats SimulateProcess(const std::string& out,
                                std::vector<double>::const_iterator from,
                                std::vector<double>::const_iterator to,
                                ctr::PerAggregateMeanEstimator* estimator) {
  SimulationStats stats;

  bool plot = !out.empty();
  std::vector<std::pair<double, double>> demand_to_plot;
  std::vector<std::pair<double, double>> leftover_demand_to_plot;

  // Ranges for which the process fits.
  std::vector<nc::viz::ColoredRange> fit_ranges;

  bool currently_fits = true;
  double no_fit_range_start = 0;

  // Measurements so far.
  std::vector<double> so_far = {*from};

  // Estimates so far.
  std::vector<double> estimates = {*from};

  double leftover_demand = 0;
  size_t timestep_count = std::distance(from, to);

  CHECK(timestep_count > 0);
  for (size_t i = 0; i < timestep_count; ++i) {
    double timestep_demand = *std::next(from, i);
    stats.total_demand += timestep_demand;

    if (plot) {
      demand_to_plot.emplace_back(i, timestep_demand);
      leftover_demand_to_plot.emplace_back(i, leftover_demand);
    }

    leftover_demand += timestep_demand;
    double estimate = estimator->MeanForNextTimestep(so_far, estimates);
    estimates.emplace_back(estimate);

    if (timestep_demand > estimate) {
      stats.total_unmet_demand += timestep_demand - estimate;
    }

    if (estimate < leftover_demand) {
      ++stats.timesteps_with_unment_demand;
      so_far.emplace_back(estimate);

      if (plot && currently_fits) {
        currently_fits = false;
        no_fit_range_start = i;
      }
    } else {
      so_far.emplace_back(leftover_demand);
      stats.wasted_capacity += estimate - leftover_demand;
      if (plot && !currently_fits) {
        currently_fits = true;
        fit_ranges.push_back({no_fit_range_start, static_cast<double>(i)});
      }
    }

    leftover_demand -= so_far.back();
    if (!FLAGS_cumulative) {
      leftover_demand = 0;
    }
  }

  stats.next_value_estimate = estimator->MeanForNextTimestep(so_far, estimates);

  if (plot) {
    if (!currently_fits) {
      fit_ranges.push_back(
          {no_fit_range_start, static_cast<double>(timestep_count)});
    }

    std::vector<std::pair<double, double>> estimate_to_plot;
    std::vector<std::pair<double, double>> signal_to_plot;
    for (size_t i = 0; i < so_far.size(); ++i) {
      signal_to_plot.push_back({i + 1, so_far[i]});
      estimate_to_plot.push_back({i, estimates[i]});
    }

    nc::viz::PlotParameters2D params;
    params.ranges = {fit_ranges};
    nc::viz::LinePlot plot(params);
    plot.AddData("Demand", demand_to_plot);
    plot.AddData("Short", leftover_demand_to_plot);
    plot.AddData("Estimate", estimate_to_plot);
    plot.AddData("Signal", signal_to_plot);
    plot.PlotToDir(out);
  }

  return stats;
}

std::vector<SimulationStats> TestBaseline(
    const std::vector<std::vector<double>>& time_series) {
  std::vector<SimulationStats> out;
  for (uint32_t i = 0; i < time_series.size(); ++i) {
    std::string plot_output = nc::StrCat("baseline_plot_", i);
    ctr::BaselineEstimator baseline_est(time_series[i]);
    SimulationStats stats =
        SimulateProcess(plot_output, time_series[i].begin(),
                        time_series[i].end(), &baseline_est);
    out.emplace_back(stats);
  }

  return out;
}

std::vector<std::vector<double>> TimeSeriesFromFile(const std::string& file) {
  std::vector<std::vector<double>> out;
  nc::File::ReadLines(file, [&out](const std::string& line) {
    std::vector<double> series;
    for (const std::string& level_str : nc::Split(line, ",")) {
      double level;
      CHECK(nc::safe_strtod(level_str, &level));
      series.emplace_back(level);
    }

    out.emplace_back(series);
  });

  return out;
}

std::vector<SimulationStats> TestEstimator(
    const std::vector<std::vector<double>>& time_series) {
  std::vector<SimulationStats> out;
  for (uint32_t i = 0; i < time_series.size(); ++i) {
    std::string plot_output = nc::StrCat("ts_plot_", i);
    const std::vector<double>& ts = time_series[i];

    double max = 1;
    if (FLAGS_simple) {
      max = FLAGS_decay_fraction;
    }

    ctr::MeanScaleEstimatorConfig cfg(FLAGS_fixed_scale, FLAGS_decay_fraction,
                                      max, FLAGS_window_size);
    ctr::MeanScaleEstimator mse(cfg);
    SimulationStats stats =
        SimulateProcess(plot_output, ts.begin(), ts.end(), &mse);
    out.emplace_back(stats);
  }

  return out;
}

static void TestConfig(const std::vector<std::vector<double>>& time_series) {
  std::vector<SimulationStats> out = TestEstimator(time_series);
  for (SimulationStats stats : out) {
    LOG(INFO) << stats.ToString();
  }
}

int main(int argc, char** argv) {
  using namespace std::chrono;
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_time_series_file.empty());

  std::vector<std::vector<double>> series =
      TimeSeriesFromFile(FLAGS_time_series_file);
  std::reverse(series.begin(), series.end());
  //  series = {series[3]};
  TestConfig(series);

  LOG(ERROR) << "Baseline:";
  std::vector<SimulationStats> out = TestBaseline(series);
  for (SimulationStats stats : out) {
    LOG(INFO) << stats.ToString();
  }
}
