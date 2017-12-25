#include <gflags/gflags.h>
#include <iostream>
#include <string>
#include <vector>

#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"

DEFINE_string(topology_root, "", "Root for topologies. Required.");
DEFINE_double(link_delay_scale, 1.0,
              "All link delays will be scaled by this number");
DEFINE_double(link_capacity_scale, 1.0,
              "All link bandwidths will be scaled by this number");

static constexpr char kTopologyExtension[] = ".graph";

static std::vector<std::string> GetTopologyFiles() {
  std::string top_root = FLAGS_topology_root;
  CHECK(top_root != "");
  return nc::File::GetFilesWithExtension(top_root, kTopologyExtension);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<std::string> topology_files = GetTopologyFiles();
  for (const std::string& topology_file : topology_files) {
    std::vector<std::string> node_order;
    nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(topology_file), &node_order);
    builder.ScaleCapacity(FLAGS_link_capacity_scale);
    builder.ScaleDelay(FLAGS_link_delay_scale);

    std::string serialized = builder.ToRepetita(node_order);
    nc::File::WriteStringToFileOrDie(serialized, topology_file);
    LOG(INFO) << "Overwrote " << topology_file << " capacity scale "
              << FLAGS_link_capacity_scale << " delay scale "
              << FLAGS_link_delay_scale;
  }
}
