#ifndef CTR_NET_MOCK_H
#define CTR_NET_MOCK_H

#include <stddef.h>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/packet.h"
#include "ncode_common/src/logging.h"
#include "pcap_data.pb.h"

namespace ctr {

class NetMock {
  // Generates the input to the system.
  std::map<AggregateId, AggregateHistory> GenerateInput(
      const std::map<AggregateId, BinSequence>& period_sequences) const {
    std::map<AggregateId, AggregateHistory> input;
    for (const auto& aggregate_and_bins : period_sequences) {
      const AggregateId& aggregate = aggregate_and_bins.first;
      const BinSequence& bins = aggregate_and_bins.second;
      input[aggregate] = bins.GenerateHistory(history_bin_size_);
    }

    return input;
  }

  // For each link returns what the residuals will be for an optimization
  // output.
  nc::net::GraphLinkMap<std::vector<double>> CheckOutput(
      const std::map<AggregateId, BinSequence>& period_sequences,
      const RoutingConfiguration& configuration) const {
    nc::net::GraphLinkMap<BinSequence> link_to_bins;

    // First need to figure out which paths cross each link. Will also build a
    // map from paths to aggregates and path indices.
    for (const auto& aggregate_and_routes : configuration.routes()) {
      const AggregateId& aggregate = aggregate_and_routes.first;
      const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;

      std::vector<double> fractions;
      for (const auto& route_and_fraction : routes) {
        fractions.emplace_back(route_and_fraction.second);
      }

      // For each of the aggregate's paths, the bins that go on that path.
      std::vector<BinSequence> aggregate_split =
          nc::FindOrDieNoPrint(period_sequences, aggregate)
              .SplitOrDie(fractions);

      for (size_t i = 0; i < routes.size(); ++i) {
        const nc::net::Walk* path = routes[i].first;
        for (nc::net::GraphLinkIndex link : path->links()) {
          link_to_bins[link].Combine(aggregate_split[i]);
        }
      }
    }

    nc::net::GraphLinkMap<std::vector<double>> out;
    for (const auto& link_and_bins : link_to_bins) {
      nc::net::GraphLinkIndex link = link_and_bins.first;
      nc::net::Bandwidth rate = graph_->GetLink(link)->bandwidth();
      out[link].emplace_back(link_and_bins.second->Residuals(rate));
    }

    return out;
  }

  std::map<AggregateId, BinSequence> GetNthPeriod(size_t n) const {
    size_t period_start_bin = n * period_duration_bins_;
    size_t period_end_bin = (n + 1) * period_duration_bins_;

    std::map<AggregateId, BinSequence> out;
    for (const auto& aggregate_and_bins : initial_sequences_) {
      const AggregateId& aggregate = aggregate_and_bins.first;
      const BinSequence& bins = aggregate_and_bins.second;
      out[aggregate] = bins.LimitRange(period_start_bin, period_end_bin);
    }

    return out;
  }

  std::unique_ptr<RoutingConfiguration> InitialOutput() const {
    std::map<AggregateId, AggregateHistory> input =
        GenerateInput(GetNthPeriod(0));
    return routing_system_->Update(input);
  }

  void Run() {
    std::unique_ptr<RoutingConfiguration> output = InitialOutput();
    for (size_t i = 0; i < period_count_; ++i) {
      std::map<AggregateId, BinSequence> period_sequences = GetNthPeriod(i);
      nc::net::GraphLinkMap<std::vector<double>> per_link_residuals =
          CheckOutput(period_sequences, *output);
      for (auto& link_and_residuals : per_link_residuals) {
        nc::net::GraphLinkIndex link = link_and_residuals.first;
        std::vector<double>& residuals = *link_and_residuals.second;

        std::vector<double>& all_residuals = all_residuals_[link];
        all_residuals.insert(all_residuals.end(), residuals.begin(),
                             residuals.end());
      }

      std::map<AggregateId, AggregateHistory> input =
          GenerateInput(period_sequences);
      output = routing_system_->Update(input);
    }
  }

  size_t period_count_;

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
