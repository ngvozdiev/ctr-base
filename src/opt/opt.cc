#include "opt.h"

#include <stddef.h>
#include <initializer_list>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/lp/lp.h"
#include "ncode_common/src/lp/mc_flow.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/algorithm.h"
#include "ncode_common/src/perfect_hash.h"

namespace ctr {

std::unique_ptr<RoutingConfiguration> ShortestPathOptimizer::Optimize(
    const TrafficMatrix& tm) {
  auto out = nc::make_unique<RoutingConfiguration>(tm);
  for (const auto& aggregate_and_demand : tm.demands()) {
    const AggregateId& aggregate_id = aggregate_and_demand.first;

    const nc::net::Walk* path =
        path_provider_->AvoidingPathOrNull(aggregate_id, {});
    CHECK(path != nullptr) << "Unable to find a path for aggregate";
    out->AddRouteAndFraction(aggregate_id, {{path, 1.0}});
  }

  return out;
}

static nc::net::GraphLinkMap<double> GetCapacities(
    const nc::net::GraphStorage& graph, double multiplier) {
  nc::net::GraphLinkMap<double> out;
  for (nc::net::GraphLinkIndex link : graph.AllLinks()) {
    double capacity = graph.GetLink(link)->bandwidth().Mbps();
    out[link] = capacity * multiplier;
  }

  return out;
}

class MinMaxProblem : public nc::lp::SingleCommodityFlowProblem {
 public:
  MinMaxProblem(const nc::net::GraphStorage* graph, double capacity_multiplier,
                bool also_minimize_delay)
      : nc::lp::SingleCommodityFlowProblem(
            GetCapacities(*graph, capacity_multiplier), graph),
        also_minimize_delay_(also_minimize_delay),
        capacity_multiplier_(capacity_multiplier) {}

  double Solve(
      std::map<nc::lp::SrcAndDst, std::vector<nc::lp::FlowAndPath>>* paths) {
    using namespace nc::lp;

    Problem problem(MINIMIZE);
    std::vector<ProblemMatrixElement> problem_matrix;
    VarMap link_to_variables =
        GetLinkToVariableMap(&problem, &problem_matrix, false);
    AddFlowConservationConstraints(link_to_variables, &problem,
                                   &problem_matrix);

    // There will be one max utilization variable.
    VariableIndex max_utilization_var = problem.AddVariable();
    problem.SetVariableRange(max_utilization_var, 0, Problem::kInifinity);
    problem.SetObjectiveCoefficient(max_utilization_var, 100000.0);

    // Will add per-link constraints to make each link's utilization less than
    // the max utilization.
    for (const auto& link_and_variables : link_to_variables) {
      nc::net::GraphLinkIndex link_index = link_and_variables.first;
      const nc::net::GraphLink* link = graph_->GetLink(link_index);
      const nc::net::GraphNodeMap<VariableIndex>& variables =
          *link_and_variables.second;

      // All utilization variables will be less than max utilization.
      ConstraintIndex utilization_var_constraint = problem.AddConstraint();
      problem.SetConstraintRange(utilization_var_constraint,
                                 Problem::kNegativeInifinity, 0);
      problem_matrix.emplace_back(utilization_var_constraint,
                                  max_utilization_var, -1.0);

      VariableIndex link_utilization_var = problem.AddVariable();
      problem.SetVariableRange(link_utilization_var, 0, Problem::kInifinity);
      problem_matrix.emplace_back(utilization_var_constraint,
                                  link_utilization_var, 1.0);

      ConstraintIndex constraint = problem.AddConstraint();
      problem.SetConstraintRange(constraint, 0, 0);
      problem_matrix.emplace_back(constraint, link_utilization_var, -1.0);
      double link_capacity = link->bandwidth().Mbps() * capacity_multiplier_;

      // The coefficient of each variable that goes over the link should be 1 /
      // link capacity.
      for (const auto& dst_and_variable : variables) {
        VariableIndex variable = *dst_and_variable.second;
        problem_matrix.emplace_back(constraint, variable, 1.0 / link_capacity);
      }

      double link_weight =
          std::chrono::duration<double, std::milli>(link->delay()).count();
      if (also_minimize_delay_) {
        problem.SetObjectiveCoefficient(link_utilization_var, link_weight);
      } else {
        problem.SetObjectiveCoefficient(link_utilization_var, -link_weight);
      }
    }

    // Solve the problem.
    problem.SetMatrix(problem_matrix);
    std::unique_ptr<Solution> solution = problem.Solve();
    bool solution_found = (solution->type() == nc::lp::OPTIMAL ||
                           solution->type() == nc::lp::FEASIBLE);
    if (!solution_found) {
      return std::numeric_limits<double>::max();
    }

    *paths = RecoverPathsFromSolution(link_to_variables, *solution);
    return solution->ObjectiveValue();
  }

  bool also_minimize_delay_;

  double capacity_multiplier_;
};

std::unique_ptr<RoutingConfiguration> MinMaxOptimizer::Optimize(
    const TrafficMatrix& tm) {
  MinMaxProblem problem(graph_, link_capacity_multiplier_,
                        also_minimize_delay_);

  for (const auto& aggregate_and_demand : tm.demands()) {
    const AggregateId& aggregate_id = aggregate_and_demand.first;
    nc::net::Bandwidth demand = aggregate_and_demand.second.first;
    problem.AddDemand(aggregate_id.src(), aggregate_id.dst(), demand.Mbps());
  }

  std::map<nc::lp::SrcAndDst, std::vector<nc::lp::FlowAndPath>> paths;
  problem.Solve(&paths);

  auto out = nc::make_unique<RoutingConfiguration>(tm);
  for (const auto& aggregate_and_demand : tm.demands()) {
    const AggregateId& aggregate_id = aggregate_and_demand.first;
    double demand = aggregate_and_demand.second.first.Mbps();

    // All input aggregates should have an entry in the solution.
    std::vector<nc::lp::FlowAndPath>& flow_and_paths =
        nc::FindOrDieNoPrint(paths, {aggregate_id.src(), aggregate_id.dst()});

    std::vector<RouteAndFraction> routes_for_aggregate;
    for (nc::lp::FlowAndPath& flow_and_path : flow_and_paths) {
      // Compute the fraction of the total demand that this path handles.
      double fraction = flow_and_path.flow() / demand;
      if (fraction == 0) {
        continue;
      }

      // Need to take ownership of the path and transfer it to the path
      // provider, the path provider will also return the same pointer for the
      // same path, avoiding creating new path objects for the same path every
      // time this runs.
      auto path = flow_and_path.TakeOwnershipOfPath();
      routes_for_aggregate.emplace_back(
          path_provider_->TakeOwnership(std::move(path)), fraction);
    }

    out->AddRouteAndFraction(aggregate_id, routes_for_aggregate);
  }

  return out;
}

class B4AggregateState;

// State to associate with each link in B4.
struct B4LinkState {
 public:
  B4LinkState(nc::net::GraphLinkIndex link, double capacity)
      : link_(link), remaining_capacity_(capacity) {}

  // How much fair share would it take for all non-frozen aggregates over this
  // link to congest it.
  double FairShareToCongest();

  // Advances fair share and assigns capacity to aggregates.
  void AdvanceByFairShare(double fair_share);

  void AddAggregate(B4AggregateState* aggregate_state) {
    aggregates_over_link_.insert(aggregate_state);
  }

  const std::set<B4AggregateState*>& aggregates_over_link() const {
    return aggregates_over_link_;
  }

  nc::net::GraphLinkIndex link() const { return link_; }

  bool IsCongested() const { return remaining_capacity_ < 0.1; }

 private:
  nc::net::GraphLinkIndex link_;

  // Aggregates that go over this link.
  std::set<B4AggregateState*> aggregates_over_link_;

  // Capacity left on this link.
  double remaining_capacity_;
};

// State to associate with each aggregate.
class B4AggregateState {
 public:
  B4AggregateState(const AggregateId& aggregate_id,
                   const DemandAndFlowCount& demand_and_flow_count,
                   double fair_share)
      : aggregate_id_(aggregate_id),
        demand_(demand_and_flow_count.first),
        current_path_(nullptr),
        fair_share_(fair_share) {
    fair_share_ratio_ = demand_.bps() / fair_share_;
  }

  const AggregateId& aggregate_id() const { return aggregate_id_; }

  const nc::net::Walk* current_path() const { return current_path_; }
  void set_current_path(const nc::net::Walk* path) { current_path_ = path; }

  void freeze() { frozen_ = true; }
  bool frozen() const { return frozen_; }

  double fair_share_ratio() const { return fair_share_ratio_; }

  double fair_share() const { return fair_share_; }

  // Adds some capacity to the current path.
  void AdvanceByFairShare(double fair_share) {
    CHECK(current_path_ != nullptr);
    path_to_capacity_[current_path_] += fair_share_ratio_ * fair_share;
  }

  // Returns a map from path to the capacity that will be sent over it.
  const std::map<const nc::net::Walk*, double>& path_to_capacity() const {
    return path_to_capacity_;
  }

 private:
  // Identifies the aggregate.
  AggregateId aggregate_id_;

  // Total volume in the entire aggregate.
  nc::net::Bandwidth demand_;

  // All links that are on this aggregate's current path.
  const nc::net::Walk* current_path_;

  // Allocation of paths to capacity on those paths.
  std::map<const nc::net::Walk*, double> path_to_capacity_;

  // Set to true when the aggregate reaches capacity.
  bool frozen_ = false;

  // The ratio C / F, where C is the total required demand of the aggregate and
  // F is the fair share at which the aggregate achieves its total demand.
  double fair_share_ratio_;

  // This aggregate's fair share.
  double fair_share_;
};

double B4LinkState::FairShareToCongest() {
  double sum = 0;  // sum of fair share ratios
  for (const B4AggregateState* aggregate : aggregates_over_link_) {
    if (aggregate->frozen()) {
      continue;
    }

    sum += aggregate->fair_share_ratio();
  }

  return remaining_capacity_ / sum;
}

void B4LinkState::AdvanceByFairShare(double fair_share) {
  double sum = 0;  // sum of capacities.
  for (B4AggregateState* aggregate : aggregates_over_link_) {
    if (aggregate->frozen()) {
      continue;
    }

    double capacity = aggregate->fair_share_ratio() * fair_share;
    sum += capacity;
  }

  remaining_capacity_ -= sum;
}

std::unique_ptr<RoutingConfiguration> B4Optimizer::Optimize(
    const TrafficMatrix& tm) {
  std::vector<B4AggregateState> aggregate_states;
  std::map<nc::net::GraphLinkIndex, B4LinkState> link_states;

  // Populate link states.
  for (nc::net::GraphLinkIndex link_index : graph_->AllLinks()) {
    const nc::net::GraphLink* link = graph_->GetLink(link_index);

    link_states.emplace(
        std::piecewise_construct, std::forward_as_tuple(link_index),
        std::forward_as_tuple(
            link_index, link->bandwidth().bps() * link_capacity_multiplier_));
  }

  // A constraint to avoid congested links.
  nc::net::GraphLinkSet to_avoid;

  // Populate aggregate states and initial paths.
  aggregate_states.reserve(tm.demands().size());
  for (const auto& aggregate_and_demand : tm.demands()) {
    const AggregateId& aggregate_id = aggregate_and_demand.first;
    const DemandAndFlowCount& demand_and_flow_count =
        aggregate_and_demand.second;

    double fair_share = 1.0;
    if (flow_count_as_fair_share_) {
      fair_share = demand_and_flow_count.second;
    }

    aggregate_states.emplace_back(aggregate_id, demand_and_flow_count,
                                  fair_share);
    B4AggregateState& aggregate_state = aggregate_states.back();

    const nc::net::Walk* path =
        path_provider_->AvoidingPathOrNull(aggregate_id, to_avoid);
    aggregate_state.set_current_path(path);

    for (nc::net::GraphLinkIndex link : path->links()) {
      B4LinkState& link_state = nc::FindOrDie(link_states, link);
      link_state.AddAggregate(&aggregate_state);
    }
  }

  size_t satisfied_aggregates = 0;
  double current_fair_share = 0;
  while (true) {
    // Figure out at which fair share the next event is -- either a link is
    // congested or an aggregate is satisfied.
    double fair_share_of_next_event = std::numeric_limits<double>::max();
    const B4LinkState* link_to_congest = nullptr;
    std::vector<B4AggregateState*> aggregates_to_satisfy;

    for (auto& link_and_state : link_states) {
      B4LinkState* link_state = &link_and_state.second;
      if (link_state->IsCongested()) {
        continue;
      }

      double fair_share_to_congest =
          current_fair_share + link_state->FairShareToCongest();
      if (fair_share_to_congest < fair_share_of_next_event) {
        fair_share_of_next_event = fair_share_to_congest;
        link_to_congest = link_state;
      }
    }

    for (B4AggregateState& aggregate_state : aggregate_states) {
      if (aggregate_state.frozen()) {
        continue;
      }

      double fair_share_to_satisfy = aggregate_state.fair_share();
      if (fair_share_to_satisfy < fair_share_of_next_event) {
        fair_share_of_next_event = fair_share_to_satisfy;
        aggregates_to_satisfy = {&aggregate_state};
        link_to_congest = nullptr;
      } else if (fair_share_to_satisfy == fair_share_of_next_event) {
        aggregates_to_satisfy.emplace_back(&aggregate_state);
        link_to_congest = nullptr;
      }
    }

    if (fair_share_of_next_event == std::numeric_limits<double>::max()) {
      // Done -- all links have been congested and all aggregates satisfied (or
      // we ran out of paths).
      break;
    }

    double fair_share_delta = fair_share_of_next_event - current_fair_share;
    for (auto& link_and_state : link_states) {
      B4LinkState* link_state = &link_and_state.second;
      link_state->AdvanceByFairShare(fair_share_delta);
    }
    for (B4AggregateState& aggregate_state : aggregate_states) {
      if (aggregate_state.frozen()) {
        continue;
      }

      aggregate_state.AdvanceByFairShare(fair_share_delta);
    }
    current_fair_share = fair_share_of_next_event;

    if (link_to_congest != nullptr) {
      // Will advance fair share to the congesting point. All aggregates that go
      // over the link will need new paths.
      to_avoid.Insert(link_to_congest->link());

      for (B4AggregateState* aggregate_state :
           link_to_congest->aggregates_over_link()) {
        const AggregateId& aggregate_id = aggregate_state->aggregate_id();
        const nc::net::Walk* path =
            path_provider_->AvoidingPathOrNull(aggregate_id, to_avoid);
        if (path == nullptr) {
          aggregate_state->freeze();
          continue;
        }

        aggregate_state->set_current_path(path);
        for (nc::net::GraphLinkIndex link : path->links()) {
          B4LinkState& link_state = nc::FindOrDie(link_states, link);
          link_state.AddAggregate(aggregate_state);
        }
      }
    }

    for (B4AggregateState* aggregate_to_satisfy : aggregates_to_satisfy) {
      // Will freeze the aggregate, since it has satisfied its demand.
      aggregate_to_satisfy->freeze();
      ++satisfied_aggregates;
    }
  }

  auto out = nc::make_unique<RoutingConfiguration>(tm);
  for (const B4AggregateState& aggregate_state : aggregate_states) {
    const std::map<const nc::net::Walk*, double>& path_to_capacity =
        aggregate_state.path_to_capacity();

    // If we were unable to satisfy the demand the total capacity allocated to
    // paths may be less than the aggregate's total demand. Will normalize.
    double total_capacity = 0;

    for (const auto& path_and_capacity : path_to_capacity) {
      total_capacity += path_and_capacity.second;
    }

    std::vector<RouteAndFraction> routes_for_aggregate;
    for (const auto& path_and_capacity : path_to_capacity) {
      const nc::net::Walk* path = path_and_capacity.first;
      double capacity = path_and_capacity.second;
      if (capacity == 0) {
        continue;
      }

      double fraction = capacity / total_capacity;
      routes_for_aggregate.emplace_back(path, fraction);
    }
    out->AddRouteAndFraction(aggregate_state.aggregate_id(),
                             routes_for_aggregate);
  }

  return out;
}

}  // namespace ctr
