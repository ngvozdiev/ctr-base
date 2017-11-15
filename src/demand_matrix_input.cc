#include "demand_matrix_input.h"

#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <random>
#include <string>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/strutil.h"

DEFINE_string(
    tm_files, "",
    "Traffic matrix file(s). If empty will look for TMs named like the "
    "topology, but with the .demands extension");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_uint64(tm_per_topology, std::numeric_limits<uint64_t>::max(),
              "How many TMs to choose for each topology");
DEFINE_uint64(seed, 1, "Seed used when choosing TMs for each topology");

using namespace std::chrono;

namespace ctr {

static std::vector<std::string> GetMatrixFiles(
    const std::string& topology_file) {
  if (!FLAGS_tm_files.empty()) {
    std::vector<std::string> out;
    for (const std::string& tm_file : nc::Split(FLAGS_tm_files, ",")) {
      std::vector<std::string> tm_globbed = nc::Glob(tm_file);
      out.insert(out.end(), tm_globbed.begin(), tm_globbed.end());
    }

    return out;
  }

  std::string matrix_location =
      nc::StringReplace(topology_file, ".graph", ".*.demands", true);
  return nc::Glob(matrix_location);
}

std::pair<std::vector<TopologyAndFilename>,
          std::vector<DemandMatrixAndFilename>>
GetDemandMatrixInputs() {
  std::vector<TopologyAndFilename> topologies = GetTopologyInputs();

  // Inputs.
  std::vector<DemandMatrixAndFilename> demands_and_matrices;

  // Randomness to pick tms for each topology.
  std::mt19937 rnd(FLAGS_seed);

  for (const TopologyAndFilename& topology : topologies) {
    const std::vector<std::string>& node_order = topology.node_order;
    std::vector<std::string> matrix_files = GetMatrixFiles(topology.name);
    if (matrix_files.empty()) {
      LOG(ERROR) << "No matrices for " << topology;
      continue;
    }

    std::vector<DemandMatrixAndFilename> inputs_for_topology;
    for (const std::string& tm_file : matrix_files) {
      std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
          nc::lp::DemandMatrix::LoadRepetitaOrDie(
              nc::File::ReadFileToStringOrDie(tm_file), node_order,
              topology.graph.get());
      if (!demand_matrix) {
        LOG(INFO) << "Skipping " << tm_file << " inconsistent TM";
        continue;
      }

      if (demand_matrix->IsTriviallySatisfiable()) {
        LOG(INFO) << "Skipping " << tm_file << " trivially satisfiable";
        continue;
      }

      demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);
      inputs_for_topology.emplace_back(tm_file, std::move(demand_matrix));
    }

    std::shuffle(inputs_for_topology.begin(), inputs_for_topology.end(), rnd);
    inputs_for_topology.resize(
        std::min(inputs_for_topology.size(),
                 static_cast<size_t>(FLAGS_tm_per_topology)));
    for (auto& input : inputs_for_topology) {
      demands_and_matrices.emplace_back(input.demand_file,
                                        std::move(input.demand_matrix));
    }
  }

  return {std::move(topologies), std::move(demands_and_matrices)};
}

}  // namespace ctr
