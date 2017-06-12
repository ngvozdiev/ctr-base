#ifndef CTR_ROUTING_SYS_H
#define CTR_ROUTING_SYS_H

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

class RoutingSystem {
 public:
  virtual ~RoutingSystem() {}

  std::unique_ptr<RoutingConfiguration> Update(
      const std::map<AggregateId, AggregateHistory>& history) {
    // How much to scale each aggregate's mean.
    std::map<AggregateId, double> aggregate_scale;

    ProbModel model(prob_model_config_);
    for (const auto& aggregate_and_history : history) {
      const AggregateId& aggregate = aggregate_and_history.first;
      const AggregateHistory& aggregate_history = aggregate_and_history.second;

      model.AddAggregate()
    }
  }

 private:
  bool CheckWithProbModel(
      const RoutingConfiguration& routing,
      const std::map<AggregateId, AggregateHistory>& histories) {
    nc::net::GraphLinkMap<std::vector<std::pair<AggregateId, double>>>
        per_link_aggregates;
    for (const auto& aggregate_and_routes : routing.routes()) {
      const AggregateId& aggregate = aggregate_and_routes.first;
      const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;
      for (const auto& route_and_fraction : routes) {
        const nc::net::Walk* walk = route_and_fraction.first;
        double fraction = route_and_fraction.second;

        for (nc::net::GraphLinkIndex link : walk->links()) {
          per_link_aggregates[link].emplace_back(aggregate, fraction);
        }
      }
    }

    ProbModel prob_model(prob_model_config_);


  }

  ProbModelConfig prob_model_config_;

  std::unique_ptr<Optimizer> optimizer_;
};

}  // namespace ctr

#endif
