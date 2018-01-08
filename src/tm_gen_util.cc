#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "ncode/common.h"
#include "ncode/file.h"
#include "ncode/logging.h"
#include "ncode/lp/demand_matrix.h"
#include "ncode/lp/lp.h"
#include "ncode/net/algorithm.h"
#include "ncode/net/net_common.h"
#include "ncode/perfect_hash.h"
#include "ncode/strutil.h"
#include "ncode/substitute.h"
#include "ncode/thread_runner.h"
#include "topology_input.h"

DEFINE_string(output_pattern,
              "demand_matrices/scale_factor_$2/locality_$1/$0_$3.demands",
              "Traffic matrices will be saved to files named after this "
              "pattern, with $0 replaced by the topology name, $1 replaced by "
              "locality, $2 replaced by scale factor and $3 replaced by a "
              "unique integer identifier.");
DEFINE_string(min_scale_factor, "1.3",
              "The generated matrix should be scaleable by at least this much "
              "without becoming unfeasible.");
DEFINE_string(locality, "0.0",
              "How much of does locality factor into the traffic matrix.");
DEFINE_uint64(seed, 1ul, "Seed for the generated TM.");
DEFINE_uint64(tm_count, 1, "Number of TMs to generate.");
DEFINE_uint64(threads, 2, "Number of threads to use.");

using namespace std::chrono;

static void AddSumConstraints(
    const nc::net::GraphNodeMap<nc::net::Bandwidth>& demands,
    const nc::net::GraphNodeMap<std::vector<nc::lp::VariableIndex>>& vars,
    std::vector<nc::lp::ProblemMatrixElement>* problem_matrix,
    nc::lp::Problem* lp_problem) {
  for (const auto& node_and_demand : demands) {
    nc::net::GraphNodeIndex node = node_and_demand.first;
    double demand_Mbps = node_and_demand.second->Mbps();

    nc::lp::ConstraintIndex constraint = lp_problem->AddConstraint();
    lp_problem->SetConstraintRange(constraint, demand_Mbps, demand_Mbps);
    for (nc::lp::VariableIndex var : vars.GetValueOrDie(node)) {
      problem_matrix->emplace_back(constraint, var, 1.0);
    }
  }
}

// Generates a demand distribution based on either incoming traffic or incoming
// traffic and outgoing traffic levels. The distribution will be as local as
// possible (will have cost per byte minimized).
std::unique_ptr<nc::lp::DemandMatrix> LPDemandGenerate(
    const nc::lp::DemandMatrix& seed_matrix, double fraction_allowance) {
  nc::net::GraphNodeMap<nc::net::Bandwidth> demand_out;
  nc::net::GraphNodeMap<nc::net::Bandwidth> demand_in;
  nc::net::GraphNodeMap<std::vector<nc::lp::VariableIndex>> vars_in;
  nc::net::GraphNodeMap<std::vector<nc::lp::VariableIndex>> vars_out;
  std::map<std::pair<nc::net::GraphNodeIndex, nc::net::GraphNodeIndex>,
           nc::lp::VariableIndex> all_vars;

  const nc::net::GraphStorage* graph = seed_matrix.graph();
  nc::net::AllPairShortestPath sp({}, graph->AdjacencyList(), nullptr, nullptr);

  nc::lp::Problem lp_problem(nc::lp::MINIMIZE);
  for (const auto& element : seed_matrix.elements()) {
    nc::net::GraphNodeIndex src = element.src;
    nc::net::GraphNodeIndex dst = element.dst;

    demand_out[src] += element.demand;
    demand_in[dst] += element.demand;

    nc::lp::VariableIndex var = lp_problem.AddVariable();
    vars_out[src].emplace_back(var);
    vars_in[dst].emplace_back(var);

    double element_demand = element.demand.Mbps();
    double min = element_demand * (1 - fraction_allowance);
    double max = element_demand * (1 + fraction_allowance);
    min = std::max(0.0, min);
    lp_problem.SetVariableRange(var, min, max);

    double obj_c = sp.GetDistance(element.src, element.dst).count();
    lp_problem.SetObjectiveCoefficient(var, obj_c * element_demand);
    all_vars[std::make_pair(src, dst)] = var;
  }

  std::vector<nc::lp::ProblemMatrixElement> problem_matrix;
  AddSumConstraints(demand_in, vars_in, &problem_matrix, &lp_problem);
  AddSumConstraints(demand_out, vars_out, &problem_matrix, &lp_problem);
  lp_problem.SetMatrix(problem_matrix);

  std::unique_ptr<nc::lp::Solution> solution = lp_problem.Solve();
  CHECK(solution->type() == nc::lp::OPTIMAL ||
        solution->type() == nc::lp::FEASIBLE);

  std::vector<nc::lp::DemandMatrixElement> new_elements;
  for (const auto& src_and_dst_and_var : all_vars) {
    nc::net::GraphNodeIndex src = src_and_dst_and_var.first.first;
    nc::net::GraphNodeIndex dst = src_and_dst_and_var.first.second;
    double demand_Mbps = solution->VariableValue(src_and_dst_and_var.second);
    new_elements.emplace_back(
        src, dst, nc::net::Bandwidth::FromMBitsPerSecond(demand_Mbps));
  }

  return nc::make_unique<nc::lp::DemandMatrix>(new_elements, graph);
}

// Generates a DemandMatrix using a scheme based on Roughan's '93 CCR paper. The
// method there is extended to support geographic locality.
class DemandGenerator {
 public:
  DemandGenerator(const nc::net::GraphStorage* graph)
      : sp_({}, graph->AdjacencyList(), nullptr, nullptr),
        graph_(graph),
        sum_inverse_delays_squared_(0) {
    for (nc::net::GraphNodeIndex src : graph_->AllNodes()) {
      for (nc::net::GraphNodeIndex dst : graph_->AllNodes()) {
        if (src == dst) {
          continue;
        }

        double distance = sp_.GetDistance(src, dst).count();
        sum_inverse_delays_squared_ += std::pow(1.0 / distance, 2);
      }
    }
  }

  // Produces a demand matrix. The matrix is not guaranteed to the satisfiable.
  std::unique_ptr<nc::lp::DemandMatrix> Generate(nc::net::Bandwidth mean,
                                                 std::mt19937* rnd) const {
    // Will start by getting the total incoming/outgoing traffic at each node.
    // These will come from an exponential distribution with the given mean.
    std::exponential_distribution<double> dist(1.0 / mean.Mbps());
    nc::net::GraphNodeMap<double> incoming_traffic_Mbps;
    nc::net::GraphNodeMap<double> outgoing_traffic_Mbps;
    double total_Mbps = 0;
    double sum_in = 0;
    double sum_out = 0;
    for (nc::net::GraphNodeIndex node : graph_->AllNodes()) {
      double in = dist(*rnd);
      double out = dist(*rnd);
      incoming_traffic_Mbps[node] = in;
      outgoing_traffic_Mbps[node] = out;
      total_Mbps += in + out;
      sum_in += in;
      sum_out += out;
    }
    double sum_product = sum_in * sum_out;

    // Time to compute the demand for each aggregate.
    double total_p = 0;
    std::vector<nc::lp::DemandMatrixElement> elements;
    for (nc::net::GraphNodeIndex src : graph_->AllNodes()) {
      for (nc::net::GraphNodeIndex dst : graph_->AllNodes()) {
        double gravity_p = incoming_traffic_Mbps.GetValueOrDie(src) *
                           outgoing_traffic_Mbps.GetValueOrDie(dst) /
                           sum_product;
        if (src == dst) {
          total_p += gravity_p;
          continue;
        }

        double value_Mbps = gravity_p * total_Mbps;
        total_p += gravity_p;
        auto bw = nc::net::Bandwidth::FromMBitsPerSecond(value_Mbps);
        if (bw > nc::net::Bandwidth::Zero()) {
          elements.push_back({src, dst, bw});
        }
      }
    }
    CHECK(std::abs(total_p - 1.0) < 0.0001) << "Not a distribution " << total_p;

    return nc::make_unique<nc::lp::DemandMatrix>(elements, graph_);
  }

 private:
  // For the locality constraints will also need the delays of the N*(N-1)
  // shortest paths in the graph.
  nc::net::AllPairShortestPath sp_;

  // The graph.
  const nc::net::GraphStorage* graph_;

  // Sum of 1 / D_i where D_i is the delay of the shortest path of the i-th
  // IE-pair
  double sum_inverse_delays_squared_;

  DISALLOW_COPY_AND_ASSIGN(DemandGenerator);
};

using DemandMatrixVector = std::vector<std::unique_ptr<nc::lp::DemandMatrix>>;
using Input = std::pair<const ctr::TopologyAndFilename*, uint32_t>;

// Generates a single traffic matrix for a given graph.
void ProcessMatrix(const Input& input, const std::vector<double>& localities,
                   const std::vector<double>& scale_factors) {
  std::string topology_filename = input.first->file;
  topology_filename = nc::File::ExtractFileName(topology_filename);
  topology_filename = nc::Split(topology_filename, ".").front();

  const std::vector<std::string>& node_order = input.first->node_order;
  uint32_t id = input.second;

  for (double locality : localities) {
    for (double scale_factor : scale_factors) {
      std::string output_location =
          nc::Substitute(FLAGS_output_pattern.c_str(), topology_filename,
                         locality, scale_factor, id);
      if (nc::File::Exists(output_location)) {
        LOG(INFO) << "Skipping existing " << output_location;
        continue;
      }

      DemandGenerator generator(input.first->graph.get());
      std::mt19937 gen(FLAGS_seed + id);
      auto demand_matrix =
          generator.Generate(nc::net::Bandwidth::FromGBitsPerSecond(1), &gen);
      demand_matrix = LPDemandGenerate(*demand_matrix, locality);
      double csf = demand_matrix->MaxCommodityScaleFactor({}, 1.0);
      CHECK(csf != 0);
      CHECK(csf == csf);
      demand_matrix = demand_matrix->Scale(csf);

      double load = 1.0 / scale_factor;
      demand_matrix = demand_matrix->Scale(load);

      std::string directory = nc::File::ExtractDirectoryName(output_location);
      if (directory != "") {
        nc::File::RecursivelyCreateDir(directory, 0777);
      }

      LOG(INFO) << "Will write " << demand_matrix->ToString() << " to "
                << output_location;
      demand_matrix->UpdateProperty("locality",
                                    nc::ToStringMaxDecimals(locality, 2));
      demand_matrix->UpdateProperty("load", nc::ToStringMaxDecimals(load, 2));
      demand_matrix->UpdateProperty(
          "seed", nc::ToStringMaxDecimals(FLAGS_seed + id, 2));
      demand_matrix->ToRepetitaFileOrDie(node_order, output_location);
    }
  }
}

static std::vector<double> GetDoubles(const std::string& in) {
  std::vector<std::string> split = nc::Split(in, ",");
  std::vector<double> out;
  for (const std::string& piece : split) {
    double v;
    CHECK(nc::safe_strtod(piece, &v));
    out.emplace_back(v);
  }

  return out;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::vector<ctr::TopologyAndFilename> topologies = ctr::GetTopologyInputs();

  std::vector<Input> inputs;
  for (const auto& topology : topologies) {
    for (size_t i = 0; i < FLAGS_tm_count; ++i) {
      inputs.emplace_back(&topology, i);
    }
  }

  std::vector<double> localities = GetDoubles(FLAGS_locality);
  std::vector<double> scale_factors = GetDoubles(FLAGS_min_scale_factor);

  nc::RunInParallel<Input>(
      inputs, [&localities, &scale_factors](const Input& input) {
        ProcessMatrix(input, localities, scale_factors);
      }, FLAGS_threads);
}
