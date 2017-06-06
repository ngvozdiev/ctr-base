#include <gflags/gflags.h>
#include "ctr.h"

#include <stddef.h>
#include <chrono>
#include <cmath>
#include <string>
#include <utility>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/perfect_hash.h"
#include "ncode_common/src/substitute.h"
#include "ncode_common/src/lp/lp.h"
#include "path_provider.h"

DEFINE_bool(debug_ctr, false, "Print debugging information about CTR");

namespace ctr {

using nc::net::GraphLink;
using nc::lp::Problem;
using nc::lp::Solution;
using nc::lp::ProblemMatrixElement;
using nc::lp::ConstraintIndex;
using nc::lp::VariableIndex;

using PathPtr = const nc::net::Walk*;

static constexpr double kM1 = 1000000000.0;
static constexpr double kM2 = 0.001;
static constexpr double kMillion = 1000000.0;
static constexpr double kLinkFullThreshold = 1.05;

static bool HasFreeCapacity(const nc::net::GraphLinkSet& links_with_no_capacity,
                            const nc::net::Walk& path) {
  return !path.ContainsAny(links_with_no_capacity);
}

bool CTROptimizer::AddFreePaths(
    const nc::net::GraphLinkSet& links_with_no_capacity,
    const std::vector<AggregateId>& aggregates_ordered) {
  bool free_paths_added = false;

  // Number of paths that will go to the optimizer. This excludes single-path
  // aggregates.
  size_t total_path_count = 0;
  bool soft_limit_exceeded_printed = false;
  bool hard_limit_exceeded_printed = false;

  for (const AggregateId& aggregate_id : aggregates_ordered) {
    std::vector<PathPtr>& paths_in_output = path_map_[aggregate_id];
    if (paths_in_output.size()) {
      // If the longest path has free capacity over all of its
      // links, then we will skip adding additional paths to this
      // aggregate.
      if (HasFreeCapacity(links_with_no_capacity, *paths_in_output.back())) {
        continue;
      }
    }

    if (total_path_count > hard_path_limit_) {
      if (!hard_limit_exceeded_printed) {
        //        LOG(INFO) << "Hard limit " << hard_path_limit_ << "exceeded ";
        hard_limit_exceeded_printed = true;
      }

      // Will add at least the shortest path.
      if (paths_in_output.empty()) {
        PathPtr path = path_provider_->AvoidingPathOrNull(aggregate_id, {});
        CHECK(path != nullptr) << "Cannot get shortest path for aggregate";
      }

      continue;
    }

    std::vector<PathPtr> paths_up_to_free;
    CHECK(paths_in_output.size() <= per_aggregate_path_limit_);
    size_t remaining_paths = per_aggregate_path_limit_ - paths_in_output.size();
    if (total_path_count > soft_path_limit_) {
      remaining_paths = 0;
    }

    if (remaining_paths > 0) {
      size_t& start_k_index = ksp_indices_[aggregate_id];
      paths_up_to_free = path_provider_->KShortestUntilAvoidingPath(
          aggregate_id, links_with_no_capacity, start_k_index, remaining_paths);
      start_k_index += paths_up_to_free.size();
    }

    if (paths_up_to_free.empty()) {
      // Out of room at the aggregate, will directly skip to the next path that
      // avoids the links.
      PathPtr path = path_provider_->AvoidingPathOrNull(aggregate_id,
                                                        links_with_no_capacity);

      // There is no path that avoids the links, and we cannot add any more
      // paths to the aggregate. Nothing to do.
      if (path == nullptr) {
        continue;
      }

      free_paths_added = true;
      paths_in_output.emplace_back(path);
      if (FLAGS_debug_ctr) {
        CLOG(INFO, YELLOW) << "Added path (i) "
                           << path->ToStringNoPorts(*graph_);
      }
      total_path_count += paths_in_output.size();
      continue;
    }

    if (total_path_count > soft_path_limit_) {
      if (!soft_limit_exceeded_printed) {
        //        LOG(INFO) << "Soft limit " << soft_path_limit_ << " exceeded
        //        ";
        soft_limit_exceeded_printed = true;
      }

      paths_in_output.emplace_back(paths_up_to_free.back());
      if (FLAGS_debug_ctr) {
        CLOG(INFO, YELLOW) << "Added path (ii) "
                           << paths_up_to_free.back()->ToStringNoPorts(*graph_);
      }
    } else {
      paths_in_output.insert(paths_in_output.end(), paths_up_to_free.begin(),
                             paths_up_to_free.end());
      if (FLAGS_debug_ctr) {
        for (const auto& path : paths_up_to_free) {
          CLOG(INFO, YELLOW) << "Added path (iii) "
                             << path->ToStringNoPorts(*graph_);
        }
      }
    }

    total_path_count += paths_in_output.size();
    free_paths_added = true;
  }

  return free_paths_added;
}

std::unique_ptr<RoutingConfiguration> CTROptimizer::Optimize(
    const TrafficMatrix& tm) {
  std::vector<AggregateId> aggregates_ordered = PrioritizeAggregates(tm);

  auto out = nc::make_unique<RoutingConfiguration>(tm);
  OptimizePrivate(tm, aggregates_ordered, nullptr, out.get());
  return out;
}

std::unique_ptr<RoutingConfiguration> CTROptimizer::OptimizeWithPrevious(
    const TrafficMatrix& tm, const RoutingConfiguration& previous) {
  std::vector<AggregateId> aggregates_ordered = PrioritizeAggregates(tm);

  auto out = nc::make_unique<RoutingConfiguration>(tm);
  double os = OptimizePrivate(tm, aggregates_ordered, &previous, out.get());
  if (os > 1) {
    return Optimize(tm);
  }

  return out;
}

std::vector<AggregateId> CTROptimizer::PrioritizeAggregates(
    const TrafficMatrix& input) const {
  std::vector<AggregateId> all_aggregates;
  for (const auto& aggregate_and_demands : input.demands()) {
    all_aggregates.emplace_back(aggregate_and_demands.first);
  }

  std::sort(all_aggregates.begin(), all_aggregates.end(),
            [&input](const AggregateId& lhs, const AggregateId& rhs) {
              return nc::FindOrDieNoPrint(input.demands(), lhs).first >
                     nc::FindOrDieNoPrint(input.demands(), rhs).first;
            });
  return all_aggregates;
}

double CTROptimizer::OptimizePrivate(
    const TrafficMatrix& input,
    const std::vector<AggregateId>& aggregates_ordered,
    const RoutingConfiguration* base_solution, RoutingConfiguration* out) {
  nc::net::GraphLinkSet links_with_no_capacity;
  double prev_obj_value = std::numeric_limits<double>::max();
  double max_oversubscription = 0;

  std::map<AggregateId, std::vector<RouteAndFraction>> aggregate_outputs;
  size_t pass_count = 0;
  while (true) {
    if (FLAGS_debug_ctr) {
      LOG(INFO) << "New pass " << pass_count << " path map at " << &path_map_;
    }
    bool added_any_paths =
        AddFreePaths(links_with_no_capacity, aggregates_ordered);
    if (pass_count != 0 && !added_any_paths) {
      break;
    }

    CTROptimizerPass pass(&input, &path_map_, graph_, base_solution);
    RunOutput& run_output = pass.run_output();
    double obj_value = run_output.obj_value;
    CHECK(obj_value != std::numeric_limits<double>::max());

    aggregate_outputs = std::move(run_output.aggregate_outputs);
    links_with_no_capacity = pass.links_with_no_capacity();

    double delta = (obj_value - prev_obj_value) / obj_value;
    if (std::abs(delta) < 0.001) {
      break;
    }

    if (FLAGS_debug_ctr) {
      LOG(INFO) << obj_value << " vs " << prev_obj_value
                << " links with no capacity " << links_with_no_capacity.Count();
    }

    prev_obj_value = obj_value;
    max_oversubscription = run_output.max_oversubscription;
    //    if (max_oversubscription <= 1.0) {
    //      break;
    //    }

    ++pass_count;
  }

  if (FLAGS_debug_ctr) {
    CLOG(INFO, RED) << "FUBAR pass done obj " << prev_obj_value << " max os "
                    << max_oversubscription;
  }

  for (auto& aggregate_and_output : aggregate_outputs) {
    out->AddRouteAndFraction(aggregate_and_output.first,
                             aggregate_and_output.second);
  }

  return max_oversubscription;
}

struct PathAndCost {
  PathAndCost(const AggregateId& aggregate_id, VariableIndex variable,
              nc::net::Bandwidth demand, PathPtr path)
      : aggregate_id(aggregate_id),
        variable(variable),
        demand(demand),
        path(path) {}

  // The id of the aggregate.
  AggregateId aggregate_id;

  // This is the index of the variable that represents the fraction of this
  // path's aggregate total volume that is sent over this path.
  VariableIndex variable;

  // The total demand of the aggregate that this path belongs to.
  nc::net::Bandwidth demand;

  // The path.
  PathPtr path;
};

static double GetFraction(const PathAndCost& path_and_cost,
                          const Solution& solution) {
  VariableIndex path_variable = path_and_cost.variable;
  return solution.VariableValue(path_variable);
}

static double PathCostUnrounded(const nc::net::Walk* path) {
  return std::chrono::duration<double, std::milli>(path->delay()).count();
}

double CTROptimizerPass::OptimizeMinLinkOversubscription() {
  using namespace std::chrono;
  auto start_at = high_resolution_clock::now();
  size_t num_paths = 0;

  nc::net::GraphLinkMap<std::vector<PathAndCost>> link_to_paths;
  nc::net::GraphLinkMap<VariableIndex> link_to_oversubscription_variable;
  std::map<AggregateId, std::vector<PathAndCost>> aggregate_to_paths;
  CHECK(paths_->size() == input_->demands().size()) << "Inconsistent path map";

  // The main LP.
  Problem problem(nc::lp::MINIMIZE);

  // The matrix of variable coefficients.
  std::vector<ProblemMatrixElement> problem_matrix;

  size_t total_paths = 0;

  // Per-aggregate constraints.
  for (const auto& aggregate_id_and_flow_count : input_->demands()) {
    const AggregateId& aggregate_id = aggregate_id_and_flow_count.first;
    const DemandAndFlowCount& demand_and_flow_count =
        aggregate_id_and_flow_count.second;

    // Skip frozen aggregates.
    if (nc::ContainsKey(frozen_aggregates_, aggregate_id)) {
      continue;
    }

    // A per-aggregate constraint to make the variables that belong to each
    // aggregate sum up to 1.
    ConstraintIndex per_aggregate_constraint = problem.AddConstraint();
    problem.SetConstraintRange(per_aggregate_constraint, 1, 1);

    double flow_count = demand_and_flow_count.second;
    const std::vector<PathPtr>& paths =
        nc::FindOrDieNoPrint(*paths_, aggregate_id);

    CHECK(!paths.empty());
    PathPtr shortest_path = paths.front();
    double total_cap = 0;
    for (PathPtr path : paths) {
      double path_cap = PathLimitFraction(aggregate_id, path);
      total_cap += path_cap;

      if (FLAGS_debug_ctr) {
        LOG(ERROR) << path->ToStringNoPorts(*graph_) << " cap " << path_cap;
      }
      if (path_cap == 0) {
        continue;
      }

      // Each path in each aggregate will have a variable associated with it.
      VariableIndex variable = problem.AddVariable();
      problem.SetVariableRange(variable, 0, path_cap);

      // The cost, delay in seconds.
      double uniqueness_weight = kM2 * flow_count * PathCostUnrounded(path) /
                                 PathCostUnrounded(shortest_path);
      problem.SetObjectiveCoefficient(
          variable, flow_count * PathCostUnrounded(path) + uniqueness_weight);

      PathAndCost path_and_cost(aggregate_id, variable,
                                demand_and_flow_count.first, path);
      for (nc::net::GraphLinkIndex link : path->links()) {
        link_to_paths[link].emplace_back(path_and_cost);
      }

      aggregate_to_paths[aggregate_id].emplace_back(path_and_cost);
      ++num_paths;
      problem_matrix.emplace_back(per_aggregate_constraint, variable, 1.0);
    }

    CHECK(total_cap >= 1);
    total_paths += paths.size();
  }

  // There will be one max oversubscription variable.
  VariableIndex max_oversubscription_var = problem.AddVariable();
  problem.SetVariableRange(max_oversubscription_var, 1.0, Problem::kInifinity);
  problem.SetObjectiveCoefficient(max_oversubscription_var, kM1);

  // The max oversubscription will be at least 1. This means
  // we have to offset the objective function back by kM1.
  problem.SetObjectiveOffset(-1.0 * (kM1));

  // Will first add per-link constraints.
  for (const auto& link_and_paths : link_to_paths) {
    nc::net::GraphLinkIndex link_index = link_and_paths.first;
    const nc::net::GraphLink* link = graph_->GetLink(link_index);

    // A constraint for the link.
    ConstraintIndex constraint = problem.AddConstraint();
    problem.SetConstraintRange(constraint, Problem::kNegativeInifinity, 0);

    // There will be a per-link oversubscription variable. Each of them will be
    // less than max_oversubscription_var.
    VariableIndex oversubscription_var = problem.AddVariable();
    problem.SetVariableRange(oversubscription_var, 1.0, Problem::kInifinity);
    link_to_oversubscription_variable[link_index] = oversubscription_var;
    problem.SetObjectiveCoefficient(oversubscription_var, 1.0);

    // A constraint for the oversubscription variable.
    ConstraintIndex oversubscription_var_constraint = problem.AddConstraint();
    problem.SetConstraintRange(oversubscription_var_constraint,
                               Problem::kNegativeInifinity, 0);
    problem_matrix.emplace_back(oversubscription_var_constraint,
                                oversubscription_var, 1);
    problem_matrix.emplace_back(oversubscription_var_constraint,
                                max_oversubscription_var, -1);
    double link_capacity_mbps = link->bandwidth().Mbps();

    // If there is frozen capacity along the link will deduct it now.
    if (frozen_capacity_.HasValue(link_index)) {
      double frozen_capacity_bps = frozen_capacity_.GetValueOrDie(link_index);
      double frozen_capacity_mbps = frozen_capacity_bps / kMillion;
      link_capacity_mbps -= frozen_capacity_mbps;
      if (link_capacity_mbps < 0) {
        link_capacity_mbps = 0;
      }
    }

    problem_matrix.emplace_back(constraint, oversubscription_var,
                                -link_capacity_mbps);
    for (const PathAndCost& path_and_cost : *link_and_paths.second) {
      VariableIndex variable = path_and_cost.variable;
      double total_volume_mbps = path_and_cost.demand.Mbps();
      problem_matrix.emplace_back(constraint, variable, total_volume_mbps);
    }
  }

  problem.SetMatrix(problem_matrix);
  // problem.DumpToFile(ncode::Substitute("dump_it_$0.lp", iteration));
  auto problem_constructed_at = high_resolution_clock::now();
  if (FLAGS_debug_ctr) {
    milliseconds problem_construction_duration =
        duration_cast<milliseconds>(problem_constructed_at - start_at);
    LOG(INFO) << nc::Substitute(
        "Problem constructed in $0ms, non-zero matrix elements $1, paths $2 ",
        duration_cast<milliseconds>(problem_construction_duration).count(),
        problem_matrix.size(), total_paths);
  }

  std::unique_ptr<Solution> solution = problem.Solve();
  bool solution_optimal = (solution->type() == nc::lp::OPTIMAL);
  if (!solution_optimal) {
    // The solver produced and unfeasible solution.
    if (solution->type() != nc::lp::FEASIBLE) {
      return std::numeric_limits<double>::max();
    }
  }

  // If the network is oversubscibed there will be one or more links with max
  // oversubscription. There will also be at least one aggregate that will have
  // all of its paths go through the maximally oversubscribed links. Will find
  // these aggregates and freeze them.
  latest_run_max_oversubscription_ =
      solution->VariableValue(max_oversubscription_var);

  auto problem_solved_at = high_resolution_clock::now();
  if (FLAGS_debug_ctr) {
    milliseconds problem_solution_duration =
        duration_cast<milliseconds>(problem_solved_at - problem_constructed_at);
    LOG(INFO) << nc::Substitute(
        "Problem solved in $0ms obj $1 paths $2, max os $3",
        duration_cast<milliseconds>(problem_solution_duration).count(),
        solution->ObjectiveValue(), total_paths,
        latest_run_max_oversubscription_);
  }

  //  std::map<AggregateId, std::set<PathPtr>>
  //      aggregate_to_max_oversubscribed_paths;
  std::set<AggregateId> max_oversubscribed_aggregates;
  for (const auto& link_and_paths : link_to_paths) {
    nc::net::GraphLinkIndex link_index = link_and_paths.first;
    const nc::net::GraphLink* link = graph_->GetLink(link_index);

    VariableIndex oversubscription_variable =
        link_to_oversubscription_variable[link_index];
    double oversubscription =
        solution->VariableValue(oversubscription_variable);
    bool eq =
        std::fabs(oversubscription - latest_run_max_oversubscription_) < 0.001;

    nc::net::Bandwidth total_load = nc::net::Bandwidth::Zero();
    for (const PathAndCost& path_over_link : *link_and_paths.second) {
      // Want to skip paths that have no flows on them.
      double fraction = GetFraction(path_over_link, *solution);
      if (fraction == 0) {
        continue;
      }

      AggregateId aggregate_id = path_over_link.aggregate_id;
      total_load += path_over_link.demand * fraction;

      if (eq) {
        //        aggregate_to_max_oversubscribed_paths[aggregate_id].insert(
        //            path_over_link.path);
        max_oversubscribed_aggregates.emplace(aggregate_id);
      }
    }

    double load_fraction = total_load / link->bandwidth();
    if (oversubscription > 1.0 || load_fraction >= kLinkFullThreshold) {
      links_with_no_capacity_.Insert(link_index);
      if (FLAGS_debug_ctr) {
        CLOG(ERROR, GREEN) << "Link with no capacity " << link->ToString();
      }
    }
  }

  for (const AggregateId& aggregate : max_oversubscribed_aggregates) {
    if (FLAGS_debug_ctr) {
      CLOG(INFO, GREEN) << "Frozen aggregate " << aggregate.ToString(*graph_);
    }
    frozen_aggregates_.insert(aggregate);
  }

  //  for (const auto& aggregate_and_paths : aggregate_to_paths) {
  //    const AggregateId& aggregate_id = aggregate_and_paths.first;
  //    const std::vector<PathAndCost>& paths = aggregate_and_paths.second;
  //
  //    std::set<PathPtr> all_paths_in_aggregate;
  //    for (const PathAndCost& path_and_cost : paths) {
  //      if (GetFraction(path_and_cost, *solution) == 0) {
  //        continue;
  //      }
  //
  //      all_paths_in_aggregate.insert(path_and_cost.path);
  //    }
  //
  //    const std::set<PathPtr>& oversubscribed_paths =
  //        aggregate_to_max_oversubscribed_paths[aggregate_id];
  //    if (all_paths_in_aggregate == oversubscribed_paths) {
  //      if (FLAGS_debug_ctr) {
  //        CLOG(INFO, GREEN) << "Frozen aggregate "
  //                          << aggregate_id.ToString(*graph_);
  //      }
  //
  //      frozen_aggregates_.insert(aggregate_id);
  //    }
  //  }

  // Will also freeze all aggregates whose all paths go over links with no
  // capacity, as this would render the next iteration unfeasible. Have to be
  // careful to also handle aggregates that have some paths that go over links
  // with capacity, but are capped.
  for (const auto& aggregate_and_paths : aggregate_to_paths) {
    const AggregateId& aggregate_id = aggregate_and_paths.first;
    const std::vector<PathAndCost>& paths = aggregate_and_paths.second;

    bool all_bad = true;
    double total_cap = 0.0;
    for (const PathAndCost& path_and_cost : paths) {
      PathPtr path = path_and_cost.path;
      if (!path->ContainsAny(links_with_no_capacity_)) {
        all_bad = false;
        total_cap += PathLimitFraction(aggregate_id, path);
      }
    }

    if (all_bad) {
      if (FLAGS_debug_ctr) {
        CLOG(INFO, GREEN) << "Frozen aggregate "
                          << aggregate_id.ToString(*graph_)
                          << " (all paths bad)";
      }
      frozen_aggregates_.insert(aggregate_id);
      continue;
    }

    if (total_cap < 1.0) {
      if (FLAGS_debug_ctr) {
        CLOG(INFO, GREEN) << "Frozen aggregate "
                          << aggregate_id.ToString(*graph_) << " (cap too low)";
      }
      frozen_aggregates_.insert(aggregate_id);
    }
  }

  if (FLAGS_debug_ctr) {
    CLOG(INFO, GREEN) << "Frozen " << frozen_aggregates_.size() << "/"
                      << paths_->size() << " aggregates";
  }

  // Populate the outputs for all frozen aggregates.
  for (const auto& aggregate_and_paths : aggregate_to_paths) {
    const AggregateId& aggregate_id = aggregate_and_paths.first;
    const std::vector<PathAndCost>& paths = aggregate_and_paths.second;
    if (!nc::ContainsKey(frozen_aggregates_, aggregate_id)) {
      continue;
    }

    // Simple sanity check.
    double total = 0;
    for (const PathAndCost& path_and_cost : paths) {
      total += GetFraction(path_and_cost, *solution);
    }
    CHECK(std::fabs(total - 1.0) < 0.0001);

    for (const PathAndCost& path_and_cost : paths) {
      double fraction = GetFraction(path_and_cost, *solution);
      fraction = std::max(0.0, fraction);
      fraction = std::min(1.0, fraction);
      if (fraction == 0) {
        continue;
      }

      double total_volume = path_and_cost.demand.bps();
      PathPtr path = path_and_cost.path;
      double capacity_taken = (total_volume * fraction);

      for (nc::net::GraphLinkIndex link : path->links()) {
        frozen_capacity_[link] += capacity_taken;
      }

      std::vector<RouteAndFraction>& aggregate_output =
          run_output_.aggregate_outputs[aggregate_id];
      aggregate_output.emplace_back(path, fraction);
    }
  }

  if (FLAGS_debug_ctr) {
    auto solution_obtained_at = high_resolution_clock::now();
    milliseconds solution_obtain_duration =
        duration_cast<milliseconds>(solution_obtained_at - problem_solved_at);
    LOG(INFO) << nc::Substitute(
        "Solution obtained in $0ms",
        duration_cast<milliseconds>(solution_obtain_duration).count());
  }

  return solution->ObjectiveValue() - link_to_paths.Count();
}

CTROptimizerPass::CTROptimizerPass(const TrafficMatrix* input,
                                   const CTRPathMap* paths,
                                   const nc::net::GraphStorage* graph,
                                   const RoutingConfiguration* base_solution)
    : input_(input),
      paths_(paths),
      graph_(graph),
      latest_run_max_oversubscription_(0),
      base_solution_(base_solution) {
  if (base_solution == nullptr) {
    FreezeSinglePathAggregates();
  } else {
    initial_obj_ = 0;
    initial_oversubscription_ = 0;
  }

  Optimize();
}

double CTROptimizerPass::PathLimitFraction(const AggregateId& aggregate,
                                           const nc::net::Walk* path) const {
  if (base_solution_ == nullptr) {
    // The path is not limited.
    return 1.0;
  }

  const std::vector<RouteAndFraction>& previous_paths =
      base_solution_->FindRoutesOrDie(aggregate);
  CHECK(!previous_paths.empty());

  for (size_t i = 0; i < previous_paths.size() - 1; ++i) {
    // If the path is not the aggregate's longest path we will limit it to its
    // current value.
    const RouteAndFraction& route_and_fraction = previous_paths[i];
    if (route_and_fraction.first == path ||
        *route_and_fraction.first == *path) {
      return route_and_fraction.second;
    }
  }

  nc::net::Delay last_path_delay = previous_paths.back().first->delay();
  if (path->delay() >= last_path_delay) {
    return 1.0;
  }

  return 0.0;
}

void CTROptimizerPass::Optimize() {
  size_t num_passes = 0;
  while (true) {
    size_t num_frozen_before = frozen_aggregates_.size();
    double obj_value = OptimizeMinLinkOversubscription();
    CHECK(obj_value != std::numeric_limits<double>::max());

    if (num_passes == 0) {
      //      run_output_.obj_value = obj_value;
      run_output_.obj_value = std::max(obj_value, initial_obj_);
      run_output_.max_oversubscription =
          std::max(initial_oversubscription_, latest_run_max_oversubscription_);
    }

    ++num_passes;
    if (latest_run_max_oversubscription_ <= 1.0 ||
        frozen_aggregates_.size() == input_->demands().size()) {
      break;
    }

    if (frozen_aggregates_.size() <= num_frozen_before) {
      LOG(ERROR) << "Unable to freeze more aggregates";
      run_output_.obj_value = std::numeric_limits<double>::max();
      break;
    }
  }
}

void CTROptimizerPass::FreezeSinglePathAggregates() {
  size_t num_frozen = 0;
  double total_cost = 0;
  for (const auto& aggregate_and_paths : *paths_) {
    const AggregateId& aggregate = aggregate_and_paths.first;
    const std::vector<PathPtr>& paths = aggregate_and_paths.second;
    CHECK(!paths.empty());
    if (paths.size() > 1) {
      continue;
    }

    const DemandAndFlowCount& input =
        nc::FindOrDieNoPrint(input_->demands(), aggregate);
    PathPtr path = paths.front();
    double flow_count = input.second;
    total_cost += flow_count * PathCostUnrounded(path) + flow_count * kM2;

    std::vector<RouteAndFraction>& aggregate_paths =
        run_output_.aggregate_outputs[aggregate];
    aggregate_paths.emplace_back(path, 1.0);

    frozen_aggregates_.insert(aggregate);
    ++num_frozen;
    for (nc::net::GraphLinkIndex link : path->links()) {
      frozen_capacity_[link] += input.first.bps();
    }
  }

  double max_oversub = 0.0;
  for (const auto& link_index_and_frozen_capacity : frozen_capacity_) {
    nc::net::GraphLinkIndex link_index = link_index_and_frozen_capacity.first;
    double frozen_capacity = *link_index_and_frozen_capacity.second;
    const nc::net::GraphLink* link = graph_->GetLink(link_index);

    double raw_capacity = link->bandwidth().bps();
    double load_fraction = frozen_capacity / raw_capacity;
    max_oversub = std::max(max_oversub, load_fraction);

    if (load_fraction >= kLinkFullThreshold) {
      links_with_no_capacity_.Insert(link_index);

      if (FLAGS_debug_ctr) {
        CLOG(ERROR, GREEN) << "Link with no capacity (i) " << link->ToString();
      }
    }
  }
  initial_oversubscription_ = std::max(1.0, max_oversub);
  initial_obj_ = total_cost + (initial_oversubscription_ - 1) * kM1;

  if (FLAGS_debug_ctr) {
    LOG(INFO) << "Frozen " << num_frozen
              << " aggregates with 1 path, initial oversubscription "
              << initial_oversubscription_ << " initial obj " << initial_obj_;
  }
}

nc::net::Bandwidth CTRQuickOptimizer::MinFreeCapacity(
    const nc::net::Links& links,
    const nc::net::GraphLinkMap<nc::net::Bandwidth>& extra_capacities,
    const nc::net::GraphLinkMap<nc::net::Bandwidth>& slack_capacities) const {
  CHECK(!links.empty());
  nc::net::Bandwidth min_free = nc::net::Bandwidth::Max();
  for (const nc::net::GraphLinkIndex link : links) {
    nc::net::Bandwidth free =
        FreeCapacityOnLink(link, extra_capacities, slack_capacities);
    min_free = std::min(min_free, free);
  }

  return min_free;
}

nc::net::Bandwidth CTRQuickOptimizer::FreeCapacityOnLink(
    nc::net::GraphLinkIndex link,
    const nc::net::GraphLinkMap<nc::net::Bandwidth>& extra_capacities,
    const nc::net::GraphLinkMap<nc::net::Bandwidth>& slack_capacities) const {
  const nc::net::GraphLinkMap<double>& link_to_load = model_.link_to_load();
  double load = link_to_load.GetValueOrDie(link);
  CHECK(load <= 1 && load >= 0);
  nc::net::Bandwidth full_capacity = graph_->GetLink(link)->bandwidth();

  nc::net::Bandwidth total_load = full_capacity * load;
  if (extra_capacities.HasValue(link)) {
    total_load += extra_capacities.GetValueOrDie(link);
  }

  if (slack_capacities.HasValue(link)) {
    total_load -= slack_capacities.GetValueOrDie(link);
  }

  CHECK(total_load <= full_capacity);
  nc::net::Bandwidth free = full_capacity - total_load;
  return free;
}

void CTRQuickOptimizer::FindRoom(
    const nc::net::GraphLinkMap<nc::net::Bandwidth>& slack_capacities,
    nc::net::GraphLinkMap<nc::net::Bandwidth>* extra_capacities,
    std::vector<PathAndLoad>* paths, nc::net::Bandwidth* bw) const {
  for (PathAndLoad& path_and_load : *paths) {
    const nc::net::Walk* path = path_and_load.first;
    nc::net::Bandwidth free_along_path =
        MinFreeCapacity(path->links(), *extra_capacities, slack_capacities);
    if (free_along_path < nc::net::Bandwidth::FromBitsPerSecond(1)) {
      continue;
    }

    nc::net::Bandwidth to_take = std::min(*bw, free_along_path);
    for (nc::net::GraphLinkIndex link : path->links()) {
      (*extra_capacities)[link] += to_take;
    }

    path_and_load.second += to_take;
    *bw -= to_take;
    if (*bw == nc::net::Bandwidth::Zero()) {
      break;
    }
  }
}

nc::net::GraphLinkSet CTRQuickOptimizer::LinksWithNoCapacity(
    const nc::net::GraphLinkMap<nc::net::Bandwidth>& extra_capacities,
    const nc::net::GraphLinkMap<nc::net::Bandwidth>& slack_capacities) const {
  nc::net::GraphLinkSet out;
  for (const nc::net::GraphLinkIndex link : graph_->AllLinks()) {
    nc::net::Bandwidth free =
        FreeCapacityOnLink(link, extra_capacities, slack_capacities);
    if (free < nc::net::Bandwidth::FromBitsPerSecond(1)) {
      out.Insert(link);
    }
  }

  return out;
}

std::unique_ptr<RoutingConfiguration> CTRQuickOptimizer::Optimize(
    const TrafficMatrix& tm) {
  CHECK(tm.demands().size() == previous_->demands().size());

  // The model reflects the current state, will use it to get the available
  // capacity at path.
  OverSubModel model(*previous_);

  // As we add traffic to paths links will get more loaded. This extra load is
  // stored here. The opposite is true for offloading traffic---the newly freed
  // capacity is stored in 'slack_capacities'.
  nc::net::GraphLinkMap<nc::net::Bandwidth> extra_capacities;
  nc::net::GraphLinkMap<nc::net::Bandwidth> slack_capacities;

  // The output. Will update it as we go.
  auto out = nc::make_unique<RoutingConfiguration>(tm);

  // Will do this in a two-stage process. Will first free up capacity.
  for (const auto& aggregate_and_new_demands : tm.demands()) {
    const AggregateId& aggregate = aggregate_and_new_demands.first;
    nc::net::Bandwidth new_demand = aggregate_and_new_demands.second.first;
    nc::net::Bandwidth previous_demand =
        nc::FindOrDieNoPrint(previous_->demands(), aggregate).first;
    const std::vector<RouteAndFraction>& previous_routes =
        nc::FindOrDieNoPrint(previous_->routes(), aggregate);
    if (new_demand > previous_demand) {
      continue;
    }

    // The aggregate has either shrunk or remained the same. Will not do
    // anything---will simply copy over its routes to the new configuration.
    out->AddRouteAndFraction(aggregate, previous_routes);

    // Need to free up capacity for each of its paths proportional to how much
    // traffic was lost.
    for (const auto& route_and_fraciton : previous_routes) {
      double fraction = route_and_fraciton.second;
      nc::net::Bandwidth previous_path_load = previous_demand * fraction;
      nc::net::Bandwidth new_path_load = new_demand * fraction;
      for (nc::net::GraphLinkIndex link : route_and_fraciton.first->links()) {
        slack_capacities[link] += (previous_path_load - new_path_load);
      }
    }
  }

  for (const auto& aggregate_and_new_demands : tm.demands()) {
    const AggregateId& aggregate = aggregate_and_new_demands.first;
    nc::net::Bandwidth new_demand = aggregate_and_new_demands.second.first;
    nc::net::Bandwidth previous_demand =
        nc::FindOrDieNoPrint(previous_->demands(), aggregate).first;
    const std::vector<RouteAndFraction>& previous_routes =
        nc::FindOrDieNoPrint(previous_->routes(), aggregate);
    if (new_demand <= previous_demand) {
      // Handled in the loop above.
      continue;
    }

    // Need to find place for the new demand. Initially the new paths will be
    // the same as the old ones. Will convert from fraction to absolute demand
    // so that we can re-compute the fractions later.
    std::vector<PathAndLoad> new_routes;
    for (const auto& route_and_fraction : previous_routes) {
      new_routes.emplace_back(route_and_fraction.first,
                              previous_demand * route_and_fraction.second);
    }

    nc::net::Bandwidth delta = new_demand - previous_demand;
    FindRoom(slack_capacities, &extra_capacities, &new_routes, &delta);
    while (delta != nc::net::Bandwidth::Zero()) {
      // Traffic did not fit on the current paths. Will have to add extra paths.
      const nc::net::Walk* new_path = path_provider_->AvoidingPathOrNull(
          aggregate, LinksWithNoCapacity(extra_capacities, slack_capacities));
      if (new_path == nullptr) {
        // No path that avoids congestion exists.
        return {};
      }

      std::vector<PathAndLoad> extra_load = {
          {new_path, nc::net::Bandwidth::Zero()}};
      FindRoom(slack_capacities, &extra_capacities, &extra_load, &delta);

      // FindRoom should have managed to put some load on the path, since we
      // specifically chose it to avoid congestion.
      CHECK(extra_load[0].second != nc::net::Bandwidth::Zero());
      new_routes.emplace_back(extra_load[0]);
    }

    // Now new_routes contains the new paths. Will compute fractions.
    nc::net::Bandwidth total;
    for (const auto& path_and_load : new_routes) {
      total += path_and_load.second;
    }
    // Sanity-check that we managed to fit all demand.
    CHECK(std::abs(total.Mbps() - new_demand.Mbps()) < 0.001);

    std::vector<RouteAndFraction> to_update;
    for (const auto& path_and_load : new_routes) {
      double fraction = path_and_load.second / total;
      if (fraction == 0) {
        LOG(INFO) << "Ignoring zero fraction path " << path_and_load.second;
        continue;
      }

      to_update.emplace_back(path_and_load.first, fraction);
    }
    out->AddRouteAndFraction(aggregate, to_update);
  }

  return out;
}

}  // namespace nc
