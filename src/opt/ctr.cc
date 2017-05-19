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

static bool HasFreeCapacity(const nc::net::GraphLinkSet& links_with_no_capacity,
                            const nc::net::Walk& path) {
  return !path.ContainsAny(links_with_no_capacity);
}

bool CTROptimizer::AddFreePaths(
    const nc::net::GraphLinkSet& links_with_no_capacity,
    const std::vector<AggregateId>& aggregates_ordered, CTRPathMap* out) {
  bool free_paths_added = false;

  // Number of paths that will go to the optimizer. This excludes single-path
  // aggregates.
  size_t total_path_count = 0;
  bool soft_limit_exceeded_printed = false;
  bool hard_limit_exceeded_printed = false;

  for (const AggregateId& aggregate_id : aggregates_ordered) {
    std::vector<PathPtr>& paths_in_output = (*out)[aggregate_id];
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
        LOG(ERROR) << "Hard limit " << hard_path_limit_ << " exceeded";
        hard_limit_exceeded_printed = true;
      }

      // Will add at least the shortest path.
      if (paths_in_output.empty()) {
        PathPtr path = path_provider_->AvoidingPathOrNull(aggregate_id, {});
        CHECK(path != nullptr) << "Cannot get shortest path for aggregate";
      }

      continue;
    }

    bool out_of_paths = false;
    std::vector<PathPtr> paths_up_to_free;
    while (paths_in_output.size() + paths_up_to_free.size() <
           per_aggregate_path_limit_) {
      PathPtr path = path_provider_->NextShortestPathOrNull(aggregate_id);
      if (path == nullptr) {
        out_of_paths = true;
        break;
      }

      paths_up_to_free.emplace_back(path);
      if (!path->ContainsAny(links_with_no_capacity)) {
        break;
      }
    }

    if (paths_up_to_free.empty() && out_of_paths) {
      // Out of room in the aggregate, and know that there are no more shortest
      // paths. Nothing we can do.
      continue;
    }

    if (paths_up_to_free.empty() && !out_of_paths) {
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
      total_path_count += paths_in_output.size();
      continue;
    }

    if (total_path_count > soft_path_limit_) {
      if (!soft_limit_exceeded_printed) {
        LOG(ERROR) << "Soft limit " << soft_path_limit_ << " exceeded";
        soft_limit_exceeded_printed = true;
      }

      paths_in_output.emplace_back(paths_up_to_free.back());
    } else {
      paths_in_output.insert(paths_in_output.end(), paths_up_to_free.begin(),
                             paths_up_to_free.end());
    }

    total_path_count += paths_in_output.size();
    free_paths_added = true;
  }

  return free_paths_added;
}

std::unique_ptr<RoutingConfiguration> CTROptimizer::Optimize(
    const TrafficMatrix& tm) {
  std::vector<AggregateId> aggregates_ordered = PrioritizeAggregates(tm);

  auto out = nc::make_unique<RoutingConfiguration>();
  OptimizePrivate(tm, aggregates_ordered, out.get());
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

void CTROptimizer::OptimizePrivate(
    const TrafficMatrix& input,
    const std::vector<AggregateId>& aggregates_ordered,
    RoutingConfiguration* out) {
  CTRPathMap path_map;
  nc::net::GraphLinkSet links_with_no_capcity;
  double prev_obj_value = std::numeric_limits<double>::max();
  double max_oversubscription = 0;

  std::map<AggregateId, std::vector<RouteAndFraction>> aggregate_outputs;
  while (true) {
    LOG(INFO) << "New pass";
    bool added_any_paths =
        AddFreePaths(links_with_no_capcity, aggregates_ordered, &path_map);
    if (!added_any_paths) {
      break;
    }

    CTROptimizerPass pass(&input, &path_map, graph_, link_capacity_multiplier_);
    RunOutput& run_output = pass.run_output();
    double obj_value = run_output.obj_value;
    CHECK(obj_value != std::numeric_limits<double>::max());

    aggregate_outputs = std::move(run_output.aggregate_outputs);
    links_with_no_capcity = pass.links_with_no_capacity();

    double delta = (obj_value - prev_obj_value) / obj_value;
    if (obj_value > prev_obj_value) {
      CHECK(delta < 0.001) << "obj value " << obj_value << " vs "
                           << prev_obj_value << " delta " << delta;
      obj_value = prev_obj_value;
    }

    if (std::abs(delta) < 0.001) {
      break;
    }

    LOG(ERROR) << obj_value << " vs " << prev_obj_value << " lnc "
               << links_with_no_capcity.Count() << " remainder ";

    prev_obj_value = obj_value;
    max_oversubscription = run_output.max_oversubscription;
    if (max_oversubscription <= 1.0) {
      break;
    }
  }

  LOG(INFO) << "FUBAR pass done obj " << prev_obj_value << " max os "
            << max_oversubscription;

  for (auto& aggregate_and_output : aggregate_outputs) {
    out->AddRouteAndFraction(aggregate_and_output.first,
                             aggregate_and_output.second);
  }
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

static double PathCost(PathPtr path) {
  using namespace std::chrono;
  return duration_cast<duration<double>>(path->delay()).count();
}

static double GetFraction(const PathAndCost& path_and_cost,
                          const Solution& solution) {
  VariableIndex path_variable = path_and_cost.variable;
  return solution.VariableValue(path_variable);
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

    for (PathPtr path : paths) {
      // Each path in each aggregate will have a variable associated with it.
      VariableIndex variable = problem.AddVariable();
      problem.SetVariableRange(variable, 0, 1);
      //      LOG(ERROR) << "cookie " << cookie << " path " <<
      //      path->ToStringNoPorts()
      //                 << " c " << min_range << " " << max_range;

      // The cost, delay in seconds.
      double cost = PathCost(path);
      double uniqueness_weight =
          kM2 * flow_count * cost / PathCost(shortest_path);
      problem.SetObjectiveCoefficient(variable,
                                      flow_count * cost + uniqueness_weight);

      PathAndCost path_and_cost(aggregate_id, variable,
                                demand_and_flow_count.first, path);
      for (nc::net::GraphLinkIndex link : path->links()) {
        link_to_paths[link].emplace_back(path_and_cost);
      }

      aggregate_to_paths[aggregate_id].emplace_back(path_and_cost);
      ++num_paths;
      problem_matrix.emplace_back(per_aggregate_constraint, variable, 1.0);
    }

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

    double link_capacity_mbps =
        link->bandwidth().bps() / kMillion * link_capacity_multiplier_;

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
  milliseconds problem_construction_duration =
      duration_cast<milliseconds>(problem_constructed_at - start_at);
  LOG(INFO) << nc::Substitute(
      "Problem constructed in $0ms, non-zero matrix elements $1, paths $2",
      duration_cast<milliseconds>(problem_construction_duration).count(),
      problem_matrix.size(), total_paths);

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
  milliseconds problem_solution_duration =
      duration_cast<milliseconds>(problem_solved_at - problem_constructed_at);
  LOG(INFO) << nc::Substitute(
      "Problem solved in $0ms obj $1 paths $2, max os $3",
      duration_cast<milliseconds>(problem_solution_duration).count(),
      solution->ObjectiveValue(), total_paths,
      latest_run_max_oversubscription_);

  std::map<AggregateId, std::set<PathPtr>>
      aggregate_to_max_oversubscribed_paths;
  for (const auto& link_and_paths : link_to_paths) {
    nc::net::GraphLinkIndex link_index = link_and_paths.first;
    const nc::net::GraphLink* link = graph_->GetLink(link_index);

    VariableIndex oversubscription_variable =
        link_to_oversubscription_variable[link_index];
    double oversubscription =
        solution->VariableValue(oversubscription_variable);
    //    LOG(ERROR) << "L " << link->ToString() << " OS " << oversubscription;
    bool eq =
        std::fabs(oversubscription - latest_run_max_oversubscription_) < 0.001;

    uint64_t total_load = 0;
    for (const PathAndCost& path_over_link : *link_and_paths.second) {
      // Want to skip paths that have no flows on them.
      double fraction = GetFraction(path_over_link, *solution);
      if (fraction == 0) {
        continue;
      }

      AggregateId aggregate_id = path_over_link.aggregate_id;
      uint64_t total_volume = path_over_link.demand.bps();
      total_load += total_volume * fraction;

      if (eq) {
        aggregate_to_max_oversubscribed_paths[aggregate_id].insert(
            path_over_link.path);
        //        LOG(ERROR) << "    FP " <<
        //        path_over_link.path->ToStringNoPorts()
        //                   << " tv " << total_volume << " f " << fraction;
      }
    }

    double load_fraction =
        total_load / (link->bandwidth().bps() * link_capacity_multiplier_);
    if (oversubscription > 1.0001 || load_fraction >= 0.99) {
      links_with_no_capacity_.Insert(link_index);
    }
  }

  for (const auto& aggregate_and_paths : aggregate_to_paths) {
    const AggregateId& aggregate_id = aggregate_and_paths.first;
    const std::vector<PathAndCost>& paths = aggregate_and_paths.second;

    std::set<PathPtr> all_paths_in_aggregate;
    for (const PathAndCost& path_and_cost : paths) {
      if (GetFraction(path_and_cost, *solution) == 0) {
        continue;
      }

      all_paths_in_aggregate.insert(path_and_cost.path);
    }

    const std::set<PathPtr>& oversubscribed_paths =
        aggregate_to_max_oversubscribed_paths[aggregate_id];
    if (all_paths_in_aggregate == oversubscribed_paths) {
      frozen_aggregates_.insert(aggregate_id);
    }
  }

  // Will also freeze all aggregates whose all paths go over links with no
  // capacity, as this would render the next iteration unfeasible.
  for (const auto& aggregate_and_paths : aggregate_to_paths) {
    const AggregateId& aggregate_id = aggregate_and_paths.first;
    const std::vector<PathAndCost>& paths = aggregate_and_paths.second;

    bool all_bad = true;
    for (const PathAndCost& path_and_cost : paths) {
      PathPtr path = path_and_cost.path;
      if (!path->ContainsAny(links_with_no_capacity_)) {
        all_bad = false;
        break;
      }
    }

    if (all_bad) {
      frozen_aggregates_.insert(aggregate_id);
    }
  }

  LOG(INFO) << "Frozen " << frozen_aggregates_.size() << "/" << paths_->size()
            << " aggregates";

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

  auto solution_obtained_at = high_resolution_clock::now();
  milliseconds solution_obtain_duration =
      duration_cast<milliseconds>(solution_obtained_at - problem_solved_at);
  LOG(INFO) << nc::Substitute(
      "Solution obtained in $0ms",
      duration_cast<milliseconds>(solution_obtain_duration).count());

  return solution->ObjectiveValue() - link_to_paths.Count();
}

CTROptimizerPass::CTROptimizerPass(const TrafficMatrix* input,
                                   const CTRPathMap* paths,
                                   const nc::net::GraphStorage* graph,
                                   double link_capacity_multiplier)
    : input_(input),
      paths_(paths),
      graph_(graph),
      latest_run_max_oversubscription_(0),
      link_capacity_multiplier_(link_capacity_multiplier) {
  FreezeSinglePathAggregates();
  Optimize();
}

void CTROptimizerPass::Optimize() {
  size_t num_passes = 0;
  while (true) {
    size_t num_frozen_before = frozen_aggregates_.size();
    double obj_value = OptimizeMinLinkOversubscription();
    CHECK(obj_value != std::numeric_limits<double>::max());

    if (num_passes == 0) {
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
    total_cost += flow_count * PathCost(path) + flow_count * kM2;

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
    double capacity = (raw_capacity * link_capacity_multiplier_);
    double oversub = frozen_capacity / capacity;
    max_oversub = std::max(max_oversub, oversub);

    if (oversub > 1.0001) {
      links_with_no_capacity_.Insert(link_index);
    }
  }
  initial_oversubscription_ = std::max(1.0, max_oversub);
  initial_obj_ = total_cost + (initial_oversubscription_ - 1) * kM1;

  LOG(ERROR) << "Frozen " << num_frozen
             << " aggregates with 1 path, initial oversubscription "
             << initial_oversubscription_ << " initial obj " << initial_obj_;
}

}  // namespace nc
