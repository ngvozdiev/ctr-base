#ifndef CTR_OPT_H
#define CTR_OPT_H

#include "ncode_net/src/net_common.h"
#include "ncode_net/src/graph_query.h"

namespace ctr {

class Optimizer {
 public:
  Optimizer(std::unique_ptr<PathProvider> path_provider)
      : path_provider_(std::move(path_provider)) {}

  virtual ~Optimizer() {}

  virtual std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) = 0;

 protected:
  // Provides paths to the optimizer.
  std::unique_ptr<PathProvider> path_provider_;

  // The graph.
  const nc::net::GraphStorage* graph_storage_;

  // All links' capacity will be multiplied by this number.
  double link_capacity_multiplier_;
};

// Only ever assigns aggregates to their single shortest path.
class ShortestPathOptimizer : public Optimizer {
 public:
  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override {
    RoutingConfiguration out;
    for (const auto& aggregate_and_demand : tm.demands()) {
      const AggregateId& aggregate_id = aggregate_and_demand.first;

      std::vector<const nc::net::Walk*> paths =
          path_provider_->KShorestPaths(aggregate_id, 0);
      CHECK(paths.size() == 1) << "Unable to find a path for aggregate";
      out.AddRouteAndFraction(aggregate_id, {paths[0], 1.0});
    }

    return out;
  }
};

class MinMaxMCProblem : public nc::lp::MCProblem {
 public:
  MinMaxMCProblem(const nc::net::GraphStorage* graph_storage,
                  double capacity_multiplier = 1.0)
      : nc::lp::MCProblem({}, graph_storage, capacity_multiplier) {}

  double Solve(
      std::map<nc::lp::SrcAndDst, std::vector<nc::lp::FlowAndPath>>* paths) {
    using namespace nc::lp;

    Problem problem(MINIMIZE);
    std::vector<ProblemMatrixElement> problem_matrix;
    VarMap link_to_variables =
        GetLinkToVariableMap(false, &problem, &problem_matrix);
    AddFlowConservationConstraints(link_to_variables, &problem,
                                   &problem_matrix);

    // There will be one max utilization variable.
    VariableIndex max_utilization_var = problem.AddVariable();
    problem.SetVariableRange(max_utilization_var, 0, Problem::kInifinity);
    problem.SetObjectiveCoefficient(max_utilization_var, 10000.0);

    // Will add per-link constraints to make each link's utilization less than
    // the max utilization.
    for (const auto& link_and_variables : link_to_variables) {
      nc::net::GraphLinkIndex link_index = link_and_variables.first;
      const nc::net::GraphLink* link = graph_storage_->GetLink(link_index);
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
      double link_capacity = link->bandwidth().Mbps();

      // The coefficient of each variable that goes over the link should be 1 /
      // link capacity.
      for (const auto& dst_and_variable : variables) {
        VariableIndex variable = *dst_and_variable.second;
        problem_matrix.emplace_back(constraint, variable, 1.0 / link_capacity);
      }

      problem.SetObjectiveCoefficient(link_utilization_var, 1);
    }

    // Solve the problem.
    problem.SetMatrix(problem_matrix);
    std::unique_ptr<Solution> solution = problem.Solve();
    bool solution_found = (solution->type() == nc::lp::OPTIMAL ||
                           solution->type() == nc::lp::FEASIBLE);
    if (!solution_found) {
      return std::numeric_limits<double>::max();
    }

    *paths = RecoverPaths(link_to_variables, *solution);
    return solution->ObjectiveValue();
  }
};

class MinMaxOptimizer : public Optimizer {
  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override {
    MinMaxMCProblem problem(graph_storage_, link_capacity_multiplier_);

    for (const auto& aggregate_and_demand : tm.demands()) {
      const AggregateId& aggregate_id = aggregate_and_demand.first;
      nc::net::Bandwidth demand = aggregate_and_demand.second;

      nc::net::GraphNodeIndex src;
      nc::net::GraphNodeIndex dst;
      std::tie(src, dst) = aggregate_id;
      problem.AddCommodity(src, dst, demand);
    }

    std::map<nc::lp::SrcAndDst, std::vector<nc::lp::FlowAndPath>> paths;
    LOG(ERROR) << "MinMax solving MC problem";
    double min_utilization = problem.Solve(&paths);
    LOG(ERROR) << "MinMax util " << min_utilization;

    RoutingConfiguration out;
    for (const auto& aggregate_and_demand : tm.demands()) {
      const AggregateId& aggregate_id = aggregate_and_demand.first;
      nc::net::Bandwidth demand = aggregate_and_demand.second;

      // All input aggregates should have an entry in the solution.
      const nc::lp::FlowAndPath& flow_and_path =
          nc::FindOrDieNoPrint(paths, aggregate_id);

      double fraction = flow_and_path.flow() / demand;
    }

    MinMaxOutput output(duration_ms, min_utilization, -1);
    if (min_utilization != std::numeric_limits<double>::max()) {
      AggregateOutputMap output_map;

      for (const auto& cookie_and_aggregate_input : input.aggregates()) {
        uint64_t cookie = cookie_and_aggregate_input.first;
        const AggregateInput& aggregate_input =
            cookie_and_aggregate_input.second;
        const std::vector<ncode::lp::FlowAndPath>& paths_in_aggregate =
            paths[{aggregate_input.src(), aggregate_input.dst()}];

        AggregateOutput& aggregate_out =
            output_map.emplace(cookie, aggregate_input).first->second;
        for (const auto& flow_and_path : paths_in_aggregate) {
          double flow_mbps = flow_and_path.first.Mbps();
          double aggregate_input_mbps = aggregate_input.rate().Mbps();
          CHECK(flow_mbps <= aggregate_input_mbps) << flow_mbps << " vs "
                                                   << aggregate_input_mbps;
          const ncode::net::LinkSequence& link_sequence = flow_and_path.second;

          const ncode::net::GraphPath* path =
              graph_storage->PathFromLinksOrDie(link_sequence, cookie);
          double fraction = flow_mbps / aggregate_input_mbps;
          if (fraction == 0) {
            continue;
          }

          PathOutput path_output(fraction, path);
          aggregate_out.AddToPaths(path_output);
        }
      }

      output.SetAggregates(std::move(output_map));
    }
  }
};

}  // namespace ctr

#endif
