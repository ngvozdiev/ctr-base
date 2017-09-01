#ifndef CTR_OPT_H
#define CTR_OPT_H

#include <algorithm>
#include <memory>

#include "ncode_common/src/logging.h"
#include "../common.h"
#include "path_provider.h"

namespace ctr {

class Optimizer {
 public:
  explicit Optimizer(PathProvider* path_provider)
      : path_provider_(path_provider) {
    CHECK(path_provider_);
    graph_ = path_provider_->graph();
  }

  virtual ~Optimizer() {}

  virtual std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) = 0;

  const nc::net::GraphStorage* graph() { return path_provider_->graph(); }

  const PathProvider* path_provider() const { return path_provider_; }

 protected:
  // Provides paths to the optimizer.
  PathProvider* path_provider_;

  // The graph.
  const nc::net::GraphStorage* graph_;
};

// Only ever assigns aggregates to their single shortest path.
class ShortestPathOptimizer : public Optimizer {
 public:
  explicit ShortestPathOptimizer(PathProvider* path_provider)
      : Optimizer(path_provider) {}

  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override;
};

// Returns a solution that minimizes the maximum utilization of any link in the
// network. The solution will be optimal, and the paths that it uses will not
// come from the path provider (i.e., they will not conform to any policy).
// Should be considered a lower bound for how low link utilization can get with
// a given traffic matrix.
class MinMaxOptimizer : public Optimizer {
 public:
  MinMaxOptimizer(PathProvider* path_provider, double capacity_multiplier,
                  bool also_minimize_delay)
      : Optimizer(path_provider),
        link_capacity_multiplier_(capacity_multiplier),
        also_minimize_delay_(also_minimize_delay) {}

  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override;

  // All links' capacity is multiplied by this number.
  double link_capacity_multiplier_;

  // If true will minmize delay as a secondary objective. If false will maximize
  // it.
  bool also_minimize_delay_;
};

// Also computes the MinMax solution of an optimization problem, but only uses
// up to the K shortest paths.
class MinMaxPathBasedOptimizer : public Optimizer {
  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override {
    // Will build a mapping from a graph link to all paths that traverse that
    // link. For convenience each path is also paired with its aggregate.
    std::map<ncode::net::GraphLinkIndex, std::vector<MMPathAndAggregate>>
        link_to_paths;

    // The main LP.
    Problem problem(ncode::lp::MINIMIZE);

    // The matrix of variable coefficients.
    std::vector<ProblemMatrixElement> problem_matrix;
    ncode::net::GraphStorage* graph_storage = path_cache->graph_storage();

    size_t num_paths = 0;
    for (const auto& cookie_and_aggregate : input.aggregates()) {
      // Each aggregate is an IE pair.
      const AggregateInput& aggregate = cookie_and_aggregate.second;

      ncode::net::NodePairPathCache* ie_cache =
          path_cache->NodePairCache(std::make_tuple(
              aggregate.src(), aggregate.dst(), aggregate.cookie()));

      // A per-aggregate constraint to make the variables that belong to each
      // aggregate sum up to 1.
      ConstraintIndex per_aggregate_constraint = problem.AddConstraint();
      problem.SetConstraintRange(per_aggregate_constraint, 1.0, 1.0);

      std::vector<ncode::net::LinkSequence> paths =
          ie_cache->PathsKHopsFromShortest(k_hops_from_shortest_path);
      for (const ncode::net::LinkSequence& link_sequence : paths) {
        // Each path in each aggregate will have a variable associated with it.
        VariableIndex variable = problem.AddVariable();
        problem.SetVariableRange(variable, 0, Problem::kInifinity);
        problem_matrix.emplace_back(per_aggregate_constraint, variable, 1.0);

        const ncode::net::GraphPath* path = graph_storage->PathFromLinksOrDie(
            link_sequence, aggregate.cookie());

        for (ncode::net::GraphLinkIndex link : link_sequence.links()) {
          link_to_paths[link].emplace_back(variable, &aggregate, path);
        }
        ++num_paths;
      }
    }

    // There will be one max utilization variable.
    VariableIndex max_utilization_var = problem.AddVariable();
    problem.SetVariableRange(max_utilization_var, 0, 1);
    problem.SetObjectiveCoefficient(max_utilization_var, 1.0);

    // Will add per-link constraints.
    for (const auto& link_and_path : link_to_paths) {
      ncode::net::GraphLinkIndex link_index = link_and_path.first;
      const ncode::net::GraphLink* link = graph_storage->GetLink(link_index);

      // A constraint for the link.
      ConstraintIndex constraint = problem.AddConstraint();
      problem.SetConstraintRange(constraint, Problem::kNegativeInifinity, 0);

      // The max utilization variable's column is always -1.
      problem_matrix.emplace_back(constraint, max_utilization_var, -1.0);

      for (const MMPathAndAggregate& path_and_aggregate :
           link_and_path.second) {
        const AggregateInput* aggregate_input =
            path_and_aggregate.aggregate_input;
        VariableIndex variable = path_and_aggregate.variable;

        // The coefficient for the column is (total volume of the aggregate /
        // link
        // capacity).
        double link_capacity = link->bandwidth().Mbps();
        double value = aggregate_input->rate().Mbps() / link_capacity;
        problem_matrix.emplace_back(constraint, variable, value);
      }
    }

    // Solve the problem.
    problem.SetMatrix(problem_matrix);
    std::unique_ptr<Solution> solution = problem.Solve();

    bool solution_found = (solution->type() == ncode::lp::OPTIMAL ||
                           solution->type() == ncode::lp::FEASIBLE);
    if (!solution_found) {
      return std::make_pair(std::numeric_limits<double>::max(),
                            AggregateOutputMap());
    }

    double best_value = solution->ObjectiveValue();

    // Recover the solution and return it.
    std::map<uint64_t, AggregateOutput> aggregate_outputs;
    std::set<uint32_t> paths_added;
    for (const auto& link_and_path : link_to_paths) {
      for (const MMPathAndAggregate& path_and_aggregate :
           link_and_path.second) {
        const AggregateInput* aggregate_input =
            path_and_aggregate.aggregate_input;
        const GraphPath* path = path_and_aggregate.path;
        if (ncode::ContainsKey(paths_added, path->tag())) {
          continue;
        }
        paths_added.insert(path->tag());

        VariableIndex variable = path_and_aggregate.variable;
        double fraction = std::max(0.0, solution->VariableValue(variable));
        if (fraction == 0) {
          continue;
        }

        uint64_t cookie = aggregate_input->cookie();
        PathOutput path_output(fraction, path);

        auto it = aggregate_outputs.emplace(cookie, *aggregate_input).first;
        AggregateOutput& aggregate_output = it->second;
        aggregate_output.AddToPaths(path_output);
      }
    }
    return std::make_pair(best_value, aggregate_outputs);
  }
};

// Runs a heuristic similar to that of B4. Each aggregate's "fair share" will be
// set to its flow count from the TM.
class B4Optimizer : public Optimizer {
 public:
  B4Optimizer(PathProvider* path_provider, bool flow_count_as_fair_share,
              double capacity_multiplier)
      : Optimizer(path_provider),
        flow_count_as_fair_share_(flow_count_as_fair_share),
        link_capacity_multiplier_(capacity_multiplier) {}

  std::unique_ptr<RoutingConfiguration> Optimize(
      const TrafficMatrix& tm) override;

 private:
  bool flow_count_as_fair_share_;

  // All links' capacity is multiplied by this number.
  double link_capacity_multiplier_;
};

}  // namespace ctr

#endif
