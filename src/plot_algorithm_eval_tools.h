// Classes and utilities to parse the metrics generated by opt_eval_util.

#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace ctr {
namespace alg_eval {

struct AggregateTMState {
  // Delay of the shortest path.
  uint32_t sp_delay_ms;

  // Demand of the aggregate.
  double rate_Mbps;

  // For each path delay and flow count.
  std::vector<std::pair<uint32_t, uint32_t>> paths;
};

struct TMState {
  std::vector<uint32_t> aggregate_sp_delay_ms;
  std::vector<double> aggregate_rate_Mbps;
};

struct OptimizerTMState {
  std::vector<uint32_t> path_flow_count;
  std::vector<uint32_t> path_delay_ms;
  std::vector<uint32_t> aggregate_path_count;
  std::vector<double> link_utilization;
  const TMState* parent_state;

  std::vector<AggregateTMState> GetAggregates() const;
};

using TopologyAndTM = std::pair<std::string, std::string>;
using TMStateMap = std::map<TopologyAndTM, std::unique_ptr<OptimizerTMState>>;

using DataVector = std::vector<double>;
using DataMap = std::map<std::pair<std::string, std::string>, DataVector>;

// Stores data about the flows in traffic matrices.
class DataStorage {
 public:
  OptimizerTMState* GetOptimizerTMState(const std::string& topology,
                                        const std::string& optimizer,
                                        const std::string& tm);

  TMState* GetTMState(const std::string& topology, const std::string& tm);

  const std::map<std::string, TMStateMap>& data() const { return data_; }

  const std::map<TopologyAndTM, std::unique_ptr<TMState>>& tm_state_data()
      const {
    return tm_state_data_;
  }

 private:
  TMState* GetTMStatePrivate(const std::string& topology,
                             const std::string& tm);

  std::map<std::string, TMStateMap> data_;
  std::map<TopologyAndTM, std::unique_ptr<TMState>> tm_state_data_;

  // Protects 'data_'.
  std::mutex mu_;
};

// Returns the percentiles of the distribution of a per-path value and the max.
std::pair<std::vector<double>, std::vector<double>> GetStretchDistribution(
    const TMStateMap& tm_state_map);

// Same as above, but only for a single run.
std::vector<double> GetStretchDistribution(const OptimizerTMState& tm_state_map,
                                           bool absolute);

std::unique_ptr<DataStorage> ParseTMGenUtilMetrics(
    const std::string& location, const std::set<TopologyAndTM>& to_ignore = {});

}  // namespace alg_eval
}  // namespace ctr