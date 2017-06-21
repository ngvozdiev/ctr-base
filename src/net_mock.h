#ifndef CTR_NET_MOCK_H
#define CTR_NET_MOCK_H

#include <stddef.h>
#include <chrono>
#include <map>
#include <memory>
#include <vector>

#include "ncode_common/src/net/net_common.h"
#include "common.h"
#include "pcap_data.h"

namespace ctr {
class RoutingSystem;
} /* namespace ctr */

namespace ctr {

class NetMock {
 public:
  NetMock(const std::map<AggregateId, BinSequence>& initial_sequences,
          std::chrono::milliseconds period_duration,
          std::chrono::milliseconds history_bin_size,
          RoutingSystem* routing_system);

  void Run();

  const nc::net::GraphLinkMap<std::vector<double>>& all_residuals() const {
    return all_residuals_;
  }

 private:
  // Generates the input to the system.
  std::map<AggregateId, AggregateHistory> GenerateInput(
      const std::map<AggregateId, BinSequence>& period_sequences) const;

  // For each link returns what the residuals will be for an optimization
  // output.
  nc::net::GraphLinkMap<std::vector<double>> CheckOutput(
      const std::map<AggregateId, BinSequence>& period_sequences,
      const RoutingConfiguration& configuration) const;

  std::map<AggregateId, BinSequence> GetNthPeriod(size_t n) const;

  std::unique_ptr<RoutingConfiguration> InitialOutput() const;

  size_t period_count_;

  std::chrono::milliseconds period_duration_;

  size_t period_duration_bins_;

  std::chrono::milliseconds history_bin_size_;

  // For each aggregate, a series of bins.
  std::map<AggregateId, BinSequence> initial_sequences_;

  nc::net::GraphLinkMap<std::vector<double>> all_residuals_;

  RoutingSystem* routing_system_;

  const nc::net::GraphStorage* graph_;
};

}  // namespace ctr

#endif
