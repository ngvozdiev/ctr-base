#include "opt_eval.h"

#include <gflags/gflags.h>
#include <iostream>
#include <type_traits>

#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/perfect_hash.h"
#include "ncode_common/src/strutil.h"

DEFINE_string(topology_files, "", "Topology files");
DEFINE_string(
    tm_files, "",
    "Traffic matrix file(s). If empty will look for TMs named like the "
    "topology, but with the .demands extension");
DEFINE_double(link_capacity_scale, 1.0, "By how much to scale all links");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_double(delay_scale, 1.0, "By how much to scale the delays of all links");
DEFINE_uint64(topology_size_limit, 100000,
              "Topologies with size more than this will be skipped");
DEFINE_uint64(topology_delay_limit_ms, 10,
              "Topologies with diameter less than this limit will be skipped");
DEFINE_uint64(tm_per_topology, std::numeric_limits<uint64_t>::max(),
              "How many TMs to run for each topology");
DEFINE_uint64(seed, 1, "Seed used when choosing TMs for each topology");

using namespace std::chrono;

namespace ctr {

static std::vector<std::string> GetTopologyFiles() {
  std::vector<std::string> out;
  std::vector<std::string> split = nc::Split(FLAGS_topology_files, ",");
  for (const std::string& piece : split) {
    std::vector<std::string> files = nc::Glob(piece);
    out.insert(out.end(), files.begin(), files.end());
  }

  return out;
}

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

std::pair<std::vector<std::unique_ptr<nc::net::GraphStorage>>,
          std::vector<OptEvalInput>>
GetOptEvalInputs() {
  std::vector<std::string> topology_files = GetTopologyFiles();
  CHECK(!topology_files.empty());

  // All graphs.
  std::vector<std::unique_ptr<nc::net::GraphStorage>> graphs;

  // Inputs.
  std::vector<OptEvalInput> inputs;

  // Randomness to pick tms for each topology.
  std::mt19937 rnd(FLAGS_seed);

  for (const std::string& topology_file : topology_files) {
    std::vector<std::string> node_order;
    nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(topology_file), &node_order);
    builder.RemoveMultipleLinks();
    builder.ScaleCapacity(FLAGS_link_capacity_scale);
    builder.ScaleDelay(FLAGS_delay_scale);
    auto graph = nc::make_unique<nc::net::GraphStorage>(builder);

    nc::net::GraphStats stats = graph->Stats();
    CHECK(!stats.sp_delay_percentiles.empty());
    if (stats.sp_delay_percentiles.back() == nc::net::Delay::max()) {
      LOG(INFO) << "Skipping " << topology_file << " graph partitioned ";
      continue;
    }

    size_t node_count = graph->AllNodes().Count();
    if (node_count > FLAGS_topology_size_limit) {
      LOG(INFO) << "Skipping " << topology_file << " / " << topology_file
                << " size limit " << FLAGS_topology_size_limit << " vs "
                << node_count;
      continue;
    }

    nc::net::Delay diameter = graph->Stats().sp_delay_percentiles.back();
    if (diameter < milliseconds(FLAGS_topology_delay_limit_ms)) {
      LOG(INFO) << "Skipping " << topology_file << " / " << topology_file
                << " delay limit " << FLAGS_topology_delay_limit_ms
                << "ms vs diameter delay "
                << duration_cast<milliseconds>(diameter).count() << "ms";
      continue;
    }

    std::vector<std::string> matrix_files = GetMatrixFiles(topology_file);
    if (matrix_files.empty()) {
      LOG(ERROR) << "No matrices for " << topology_file;
      continue;
    }

    std::vector<OptEvalInput> inputs_for_topology;
    std::string top_file_trimmed = nc::Split(topology_file, "/").back();
    for (const std::string& tm_file : matrix_files) {
      std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
          nc::lp::DemandMatrix::LoadRepetitaOrDie(
              nc::File::ReadFileToStringOrDie(tm_file), node_order,
              graph.get());
      if (!demand_matrix) {
        LOG(INFO) << "Skipping " << tm_file << " inconsistent TM";
        continue;
      }

      demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);
      std::string tm_file_trimmed = nc::Split(tm_file, "/").back();
      inputs_for_topology.emplace_back(top_file_trimmed, tm_file_trimmed,
                                       std::move(demand_matrix));
    }

    std::shuffle(inputs_for_topology.begin(), inputs_for_topology.end(), rnd);
    inputs_for_topology.resize(
        std::min(inputs_for_topology.size(),
                 static_cast<size_t>(FLAGS_tm_per_topology)));
    for (auto& input : inputs_for_topology) {
      inputs.emplace_back(input.topology_file, input.tm_file,
                          std::move(input.demand_matrix));
    }

    graphs.emplace_back(std::move(graph));
  }

  return {std::move(graphs), std::move(inputs)};
}

}  // namespace ctr
