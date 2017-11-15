#include "topology_input.h"

#include <gflags/gflags.h>
#include <stddef.h>
#include <chrono>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/perfect_hash.h"
#include "ncode_common/src/strutil.h"

DEFINE_string(topology_files, "", "Topology files");
DEFINE_double(link_capacity_scale, 1.0, "By how much to scale all links");
DEFINE_double(delay_scale, 1.0, "By how much to scale the delays of all links");
DEFINE_uint64(topology_size_limit, 100000,
              "Topologies with size more than this will be skipped");
DEFINE_uint64(topology_delay_limit_ms, 10,
              "Topologies with diameter less than this limit will be skipped");

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

std::vector<TopologyAndFilename> GetTopologyInputs() {
  using namespace std::chrono;

  std::vector<TopologyAndFilename> out;
  std::vector<std::string> topology_files = GetTopologyFiles();
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

    out.emplace_back(std::move(graph), node_order, topology_file);
  }

  return out;
}

}  // namespace ctr
