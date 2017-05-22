#ifndef CTR_OPT_COMPARE_H
#define CTR_OPT_COMPARE_H

#include "../common.h"
#include "opt.h"

namespace ctr {

// Performance of a single optimizer on a single TM.
class PerOptimizerSummary {
 public:
  PerOptimizerSummary(const TrafficMatrix& input,
                      const RoutingConfiguration& output,
                      double link_scale_factor, uint64_t processing_time_ms);

  void SaveStats(const std::string& file_prefix) const;

  std::vector<double> LinkUtilization() const;

 private:
  // Per-flow path stretch.
  std::vector<double> path_stretches_;
  std::vector<double> path_stretches_rel_;

  // Link utilization, first is flow over link, second is link capacity.
  std::vector<std::pair<double, double>> link_loads_;

  // Number of paths per aggregate.
  std::vector<double> num_paths_;

  // Unmet demand.
  std::vector<double> unmet_demand_;

  // Time to run the optimizer.
  uint64_t processing_time_ms_;

  // Number of aggregates.
  uint64_t aggregate_count_;
};

}  // namespace ctr

#endif
