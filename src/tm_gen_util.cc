#include <gflags/gflags.h>
#include <stddef.h>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/algorithm.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/perfect_hash.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/substitute.h"
#include "ncode_common/src/thread_runner.h"
#include "opt/opt.h"
#include "topology_input.h"

DEFINE_string(output_pattern,
              "demand_matrices/scale_factor_$2/locality_$1/$0_$3.demands",
              "Traffic matrices will be saved to files named after this "
              "pattern, with $0 replaced by the topology name, $1 replaced by "
              "locality, $2 replaced by scale factor and $3 replaced by a "
              "unique integer identifier.");
DEFINE_double(min_scale_factor, 1.3,
              "The generated matrix should be scaleable by at least this much "
              "without becoming unfeasible.");
DEFINE_double(locality, 0.0,
              "How much of does locality factor into the traffic matrix.");
DEFINE_uint64(seed, 1ul, "Seed for the generated TM.");
DEFINE_uint64(tm_count, 1, "Number of TMs to generate.");
DEFINE_uint64(threads, 2, "Number of threads to use.");

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
  std::unique_ptr<nc::lp::DemandMatrix> SinglePass(double locality,
                                                   nc::net::Bandwidth mean,
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
        double gravity_p =
            (1 - locality) * incoming_traffic_Mbps.GetValueOrDie(src) *
            outgoing_traffic_Mbps.GetValueOrDie(dst) / sum_product;
        if (src == dst) {
          total_p += gravity_p;
          continue;
        }

        double distance = sp_.GetDistance(src, dst).count();
        double locality_p =
            locality / (distance * distance * sum_inverse_delays_squared_);

        double value_Mbps = (locality_p + gravity_p) * total_Mbps;
        total_p += locality_p + gravity_p;
        auto bw = nc::net::Bandwidth::FromMBitsPerSecond(value_Mbps);
        if (bw > nc::net::Bandwidth::Zero()) {
          elements.push_back({src, dst, bw});
        }
      }
    }
    CHECK(std::abs(total_p - 1.0) < 0.0001) << "Not a distribution " << total_p;

    return nc::make_unique<nc::lp::DemandMatrix>(elements, graph_);
  }

  // Returns a random matrix with the given commodity scale factor. Will
  // repeatedly call SinglePass to generate a series of matrices with the
  // highest mean rate that the commodity scale factor allows.
  std::unique_ptr<nc::lp::DemandMatrix> Generate(double commodity_scale_factor,
                                                 double locality,
                                                 std::mt19937* rnd) const {
    CHECK(commodity_scale_factor >= 1.0);
    CHECK(locality >= 0);
    CHECK(locality <= 1);

    auto demand_matrix =
        SinglePass(locality, nc::net::Bandwidth::FromMBitsPerSecond(1), rnd);
    double csf = demand_matrix->MaxCommodityScaleFactor({}, 1.0);
    CHECK(csf != 0);
    demand_matrix = demand_matrix->Scale(csf);
    demand_matrix = demand_matrix->Scale(1.0 / commodity_scale_factor);
    LOG(INFO) << "M " << csf;
    return demand_matrix;
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
void ProcessMatrix(const Input& input) {
  std::string topology_filename = input.first->file;
  topology_filename = nc::File::ExtractFileName(topology_filename);
  topology_filename = nc::Split(topology_filename, ".").front();

  const std::vector<std::string>& node_order = input.first->node_order;
  uint32_t id = input.second;

  DemandGenerator generator(input.first->graph.get());
  std::mt19937 gen(FLAGS_seed + id);
  auto demand_matrix =
      generator.Generate(FLAGS_min_scale_factor, FLAGS_locality, &gen);

  // Will ignore all aggregates less than 1Mbps.
  //  demand_matrix =
  //      demand_matrix->Filter([](const nc::lp::DemandMatrixElement& element) {
  //        return element.demand < nc::net::Bandwidth::FromMBitsPerSecond(1);
  //      });

  std::string output_location =
      nc::Substitute(FLAGS_output_pattern.c_str(), topology_filename,
                     FLAGS_locality, FLAGS_min_scale_factor, id);
  std::string directory = nc::File::ExtractDirectoryName(output_location);
  if (directory != "") {
    nc::File::RecursivelyCreateDir(directory, 0777);
  }

  LOG(INFO) << "Will write " << demand_matrix->ToString() << " to "
            << output_location;
  nc::File::WriteStringToFileOrDie(demand_matrix->ToRepetita(node_order),
                                   output_location);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::vector<ctr::TopologyAndFilename> topologies = ctr::GetTopologyInputs();

  std::vector<Input> inputs;
  uint32_t id = 0;
  for (const auto& topology : topologies) {
    for (size_t i = 0; i < FLAGS_tm_count; ++i) {
      inputs.emplace_back(&topology, id++);
    }
  }

  nc::RunInParallel<Input>(
      inputs, [](const Input& input) { ProcessMatrix(input); }, FLAGS_threads);
}
