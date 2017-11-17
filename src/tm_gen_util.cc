#include <gflags/gflags.h>
#include <stddef.h>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/substitute.h"
#include "ncode_common/src/thread_runner.h"
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

using DemandMatrixVector = std::vector<std::unique_ptr<nc::lp::DemandMatrix>>;
using Input = std::pair<const ctr::TopologyAndFilename*, uint32_t>;

// Generates a single traffic matrix for a given graph.
void ProcessMatrix(const Input& input) {
  std::string topology_filename = input.first->file;
  topology_filename = nc::File::ExtractFileName(topology_filename);
  topology_filename = nc::Split(topology_filename, ".").front();

  const std::vector<std::string>& node_order = input.first->node_order;
  uint32_t id = input.second;

  nc::lp::DemandGenerator generator(input.first->graph.get());
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
