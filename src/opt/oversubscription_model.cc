#include "oversubscription_model.h"

#include <utility>

#include "ncode/common.h"
#include "ncode/logging.h"
#include "ncode/map_util.h"

namespace ctr {

using namespace nc::net;

OverSubModelLinkState::OverSubModelLinkState(Bandwidth capacity,
                                             const GraphLink* link)
    : original_capacity_(capacity),
      link_(link),
      bottlenecked_(false),
      remaining_capacity_(capacity) {}

double OverSubModelLinkState::SubscriptionRatio() const {
  Bandwidth total_demand = Bandwidth::Zero();
  for (OverSubModelPathState* path_state : path_states_over_link_) {
    if (path_state->bottlenecked()) {
      continue;
    }

    total_demand += path_state->path_demand();
  }

  return total_demand / remaining_capacity_;
}

void OverSubModelLinkState::Bottleneck() {
  CHECK(!bottlenecked_) << "Link already bottlenecked";

  // If the link is oversubscribed all non-bottlenecked flows will get their
  // levels reduced by the oversubscription ratio.
  double ratio = SubscriptionRatio();
  if (ratio < 1.0) {
    ratio = 1.0;
  }

  for (OverSubModelPathState* path_state : path_states_over_link_) {
    if (path_state->bottlenecked()) {
      continue;
    }

    Bandwidth share = path_state->path_demand() / ratio;
    path_state->Bottleneck(share);
  }
}

void OverSubModelLinkState::AddPathState(OverSubModelPathState* path_state) {
  path_states_over_link_.emplace_back(path_state);
}

OverSubModelPathState::OverSubModelPathState(const AggregateId& aggregate_id,
                                             double fraction,
                                             nc::net::Bandwidth path_demand,
                                             const nc::net::Walk* path)
    : fraction_(fraction),
      aggregate_id_(aggregate_id),
      path_(path),
      bottlenecked_(false),
      bottleneck_rate_(Bandwidth::Zero()),
      path_demand_(path_demand) {}

void OverSubModelPathState::Bottleneck(Bandwidth at_rate) {
  CHECK(!bottlenecked_) << "Path already bottlenecked";
  if (at_rate > path_demand_) {
    at_rate = path_demand_;
  }

  bottlenecked_ = true;
  bottleneck_rate_ = at_rate;
  for (OverSubModelLinkState* link_state : link_states_in_path_) {
    link_state->ReduceCapacity(at_rate);
  }
}

Bandwidth OverSubModelPathState::bottleneck_rate() const {
  CHECK(bottlenecked_) << "Path not bottlenecked yet";
  return bottleneck_rate_;
}

void OverSubModelPathState::AddLinkState(OverSubModelLinkState* link_state) {
  link_states_in_path_.emplace_back(link_state);
}

OverSubModel::OverSubModel(const RoutingConfiguration& routing)
    : graph_(routing.graph()) {
  std::map<const nc::net::Walk*, OverSubModelPathState> path_to_state;
  std::map<nc::net::GraphLinkIndex, OverSubModelLinkState> link_to_state;
  InitState(routing, &path_to_state, &link_to_state);
  PopulateCapacities(&path_to_state, &link_to_state);

  for (const auto& path_and_state : path_to_state) {
    const Walk* path = path_and_state.first;
    const OverSubModelPathState& path_state = path_and_state.second;

    const DemandAndFlowCount& demand_and_flows =
        nc::FindOrDieNoPrint(routing.demands(), path_state.aggregate_id());
    double num_flows = demand_and_flows.second * path_state.fraction();
    per_flow_bandwidth_[path] = path_state.bottleneck_rate() / num_flows;
    if (path_state.bottleneck_rate() < path_state.path_demand()) {
      aggregates_no_fit_.insert(path_state.aggregate_id());
    }
  }

  for (const auto& link_and_state : link_to_state) {
    GraphLinkIndex link = link_and_state.first;
    const OverSubModelLinkState& link_state = link_and_state.second;
    nc::net::Bandwidth total_over_link = nc::net::Bandwidth::Zero();
    for (const OverSubModelPathState* path_state :
         link_state.path_states_over_link()) {
      total_over_link += path_state->bottleneck_rate();
    }

    link_to_load_[link] = total_over_link / link_state.original_capacity();
  }

  for (GraphLinkIndex link : graph_->AllLinks()) {
    link_to_load_[link];
  }
}

void OverSubModel::InitState(
    const RoutingConfiguration& routing,
    std::map<const nc::net::Walk*, OverSubModelPathState>* path_to_state,
    std::map<nc::net::GraphLinkIndex, OverSubModelLinkState>* link_to_state) {
  for (const auto& aggregate_and_routes : routing.routes()) {
    const AggregateId& aggregate = aggregate_and_routes.first;
    const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;

    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(routing.demands(), aggregate);
    Bandwidth demand = demand_and_flow_count.first;

    for (const auto& route_and_fraction : routes) {
      double fraction = route_and_fraction.second;
      const Walk* path = route_and_fraction.first;

      Bandwidth path_demand = demand * fraction;
      OverSubModelPathState& path_state = nc::LookupOrInsert(
          path_to_state, path,
          OverSubModelPathState(aggregate, fraction, path_demand, path));

      for (GraphLinkIndex link_index : path->links()) {
        const GraphLink* link = graph_->GetLink(link_index);

        Bandwidth capacity = link->bandwidth();
        OverSubModelLinkState& link_state = nc::LookupOrInsert(
            link_to_state, link_index, OverSubModelLinkState(capacity, link));
        path_state.AddLinkState(&link_state);
        link_state.AddPathState(&path_state);
      }
    }
  }
}

void OverSubModel::PopulateCapacities(
    std::map<const nc::net::Walk*, OverSubModelPathState>* path_to_state,
    std::map<nc::net::GraphLinkIndex, OverSubModelLinkState>* link_to_state) {
  while (true) {
    double worst_subscription_ratio = 0;
    OverSubModelLinkState* worst_subscription_ratio_link = nullptr;

    for (auto& link_and_state : *link_to_state) {
      OverSubModelLinkState& link_state = link_and_state.second;
      double subscription_ratio = link_state.SubscriptionRatio();
      if (subscription_ratio > worst_subscription_ratio) {
        worst_subscription_ratio = subscription_ratio;
        worst_subscription_ratio_link = &link_state;
      }
    }

    // If the worst subscription ratio is less than 1 no links are
    // oversubscribed -- all paths fit and we are done.
    if (worst_subscription_ratio <= 1) {
      for (auto& path_and_state : *path_to_state) {
        OverSubModelPathState& path_state = path_and_state.second;
        if (path_state.bottlenecked()) {
          continue;
        }

        path_state.Bottleneck(path_state.path_demand());
      }

      break;
    }

    // Bottleneck all paths that go over the most oversubscribed link and
    // repeat.
    worst_subscription_ratio_link->Bottleneck();
  }
}

}  // namespace fubarizer
