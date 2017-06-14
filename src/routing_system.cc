#include "routing_system.h"

#include <stddef.h>
#include <algorithm>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/net_common.h"

namespace ctr {

std::unique_ptr<RoutingConfiguration> RoutingSystem::Update(
    const std::map<AggregateId, AggregateHistory>& history) {
  // First we will predict what the history is expected to be on the next
  // timestep.
  std::map<AggregateId, AggregateHistory> next_history =
      estimator_->EstimateNext(history);

  // The initial input will assign all aggregates to their mean level.
  std::map<AggregateId, DemandAndFlowCount> input =
      GetInitialInput(next_history);

  std::unique_ptr<RoutingConfiguration> output;
  while (true) {
    TrafficMatrix tm(graph_, input);
    output = optimizer_->Optimize(tm);

    // After optimizing we should check to see which aggregates go over links
    // that do not fit.
    std::set<AggregateId> aggregates_no_fit =
        CheckWithProbModel(*output, next_history);
    if (aggregates_no_fit.empty()) {
      break;
    }

    // Need to scale up the aggregates that do not fit.
    bool scaled_any =
        ScaleUpAggregates(aggregates_no_fit, next_history, &input);
    if (!scaled_any) {
      // All aggregates at their max.
      break;
    }
  }

  CHECK(output);
  return output;
}

std::set<AggregateId> RoutingSystem::CheckWithProbModel(
    const RoutingConfiguration& routing,
    const std::map<AggregateId, AggregateHistory>& histories) {
  CHECK(histories.size() == routing.routes().size());

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

  std::vector<ProbModelQuery> queries;
  std::vector<const std::vector<std::pair<AggregateId, double>>*> map_values;
  for (const auto& link_and_aggregates : per_link_aggregates) {
    nc::net::GraphLinkIndex link = link_and_aggregates.first;
    const std::vector<std::pair<AggregateId, double>>&
        aggregates_and_fractions = *link_and_aggregates.second;
    map_values.emplace_back(&aggregates_and_fractions);

    queries.emplace_back();
    ProbModelQuery& query = queries.back();
    query.type = ProbModelQuery::BOTH;
    query.rate = graph_->GetLink(link)->bandwidth();
    query.aggregates = aggregates_and_fractions;
  }

  ProbModel prob_model(prob_model_config_);
  for (const auto& aggregate_and_history : histories) {
    const AggregateId& aggregate = aggregate_and_history.first;
    const AggregateHistory& history = aggregate_and_history.second;
    prob_model.AddAggregate(aggregate, &history);
  }
  std::vector<ProbModelReply> replies = prob_model.Query(queries);

  std::set<AggregateId> out;
  for (size_t i = 0; i < replies.size(); ++i) {
    const ProbModelReply& reply = replies[i];
    if (reply.fit) {
      continue;
    }

    for (const auto& aggregate_and_fraction : *(map_values[i])) {
      out.insert(aggregate_and_fraction.first);
    }
  }

  return out;
}

bool RoutingSystem::ScaleUpAggregates(
    const std::set<AggregateId>& aggregates,
    const std::map<AggregateId, AggregateHistory>& histories,
    std::map<AggregateId, DemandAndFlowCount>* out) {
  CHECK(!aggregates.empty());
  bool scaled_at_least_one = false;
  for (const auto& aggregate_and_history : histories) {
    const AggregateId& aggregate_id = aggregate_and_history.first;
    if (!nc::ContainsKey(aggregates, aggregate_id)) {
      continue;
    }

    const AggregateHistory& history = aggregate_and_history.second;
    nc::net::Bandwidth mean_rate = history.mean_rate();
    nc::net::Bandwidth max_rate = history.max_rate();

    const double range_mbps = max_rate.Mbps() - mean_rate.Mbps();
    if (range_mbps == 0) {
      continue;
    }

    DemandAndFlowCount& current_demand_and_flow_count =
        nc::FindOrDieNoPrint(*out, aggregate_id);
    nc::net::Bandwidth& current_rate = current_demand_and_flow_count.first;
    CHECK(current_rate >= mean_rate);

    const double mean_rate_mbps = mean_rate.Mbps();
    const double current_rate_mbps = current_rate.Mbps();
    const double curr_scale = current_rate_mbps / mean_rate_mbps;

    // Will scale the aggregate's rate to be somewhere between the mean and
    // the max. Need to first figure out where it is now.
    double curr_fraction = mean_rate_mbps * (curr_scale - 1) / range_mbps;
    if (curr_fraction > 0.99) {
      continue;
    }

    // Want to avoid updating all_fit in the case where all aggregates at
    // already at 1.0.
    scaled_at_least_one = true;

    // Will move each aggregate 10% closer to its max rate.
    double new_fraction = curr_fraction + scale_fraction_;
    new_fraction = std::min(new_fraction, 1.0);

    double new_rate_mbps = mean_rate_mbps + range_mbps * new_fraction;
    current_rate = nc::net::Bandwidth::FromMBitsPerSecond(new_rate_mbps);
  }

  return scaled_at_least_one;
}

std::map<AggregateId, DemandAndFlowCount> RoutingSystem::GetInitialInput(
    const std::map<AggregateId, AggregateHistory>& histories) {
  std::map<AggregateId, DemandAndFlowCount> out;
  for (const auto& aggregate_and_history : histories) {
    const AggregateId& aggregate_id = aggregate_and_history.first;
    const AggregateHistory& history = aggregate_and_history.second;

    out[aggregate_id] = {history.mean_rate(), history.flow_count()};
  }

  return out;
}

}  // namespace ctr
