#include <gflags/gflags.h>
#include <ncode/common.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/net/algorithm.h>
#include <ncode/net/net_common.h>
#include <ncode/net/net_gen.h>
#include <ncode/perfect_hash.h>
#include <ncode/strutil.h>
#include <stdint.h>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include "opt/ctr.h"
#include "topology_input.h"

DEFINE_double(link_delay_scale, 1.0,
              "All link delays will be scaled by this number");
DEFINE_double(link_capacity_scale, 1.0,
              "All link bandwidths will be scaled by this number");
DEFINE_uint64(split_threshold_ms, 0, "If not 0, will split the topologies");
DEFINE_bool(dump_llpd, false, "If true will dump LLP per topology and exit");
DEFINE_bool(dump_stats, false, "If true will dump overall stats and exit");
DEFINE_double(sp_fraction, 1.4, "How far from the SP a path can be");
DEFINE_double(link_fraction_limit, 0.7,
              "At least this much of the SP's links can be routed around");

static void DumpStats(const std::vector<ctr::TopologyAndFilename>& topologies) {
  std::vector<double> node_counts;
  std::vector<double> edge_counts;
  for (const auto& topology_and_filename : topologies) {
    size_t node_count = topology_and_filename.graph->NodeCount();
    size_t edge_count = topology_and_filename.graph->LinkCount();

    node_counts.emplace_back(node_count);
    edge_counts.emplace_back(edge_count / 2);
  }

  std::vector<double> np = nc::Percentiles(&node_counts, 100);
  std::vector<double> ep = nc::Percentiles(&edge_counts, 100);

  LOG(INFO) << nc::Substitute("Nodes min $0, med $1, p90: $2, max $3", np[0],
                              np[50], np[90], np[100]);
  LOG(INFO) << nc::Substitute("Edges min $0, med $1, p90: $2, max $3", ep[0],
                              ep[50], ep[90], ep[100]);
}

// Returns the links in the graph for which both the source and the destination
// nodes are in a given set of nodes.
static nc::net::GraphLinkSet LinksBetweenNodes(
    const nc::net::GraphStorage& graph, const nc::net::GraphNodeSet& nodes) {
  nc::net::GraphLinkSet out;

  const nc::net::AdjacencyList& adjacency_list = graph.AdjacencyList();
  for (nc::net::GraphNodeIndex src_node : nodes) {
    for (const nc::net::AdjacencyList::LinkInfo& link_info :
         adjacency_list.GetNeighbors(src_node)) {
      if (nodes.Contains(link_info.dst_index)) {
        out.insert(link_info.link_index);
      }
    }
  }

  return out;
}

static bool IsAlreadySeen(
    const std::vector<nc::net::GraphNodeSet>& already_seen,
    nc::net::GraphNodeIndex node) {
  for (const auto& set : already_seen) {
    if (set.Contains(node)) {
      return true;
    }
  }

  return false;
}

// Splits a graph into subgraphs by removing all links with delay larger than
// threshold. Will return a separate GraphBuilder for each connected component.
std::vector<nc::net::GraphBuilder> SplitTopology(
    const nc::net::GraphBuilder& builder, nc::net::Delay threshold) {
  nc::net::GraphStorage full_graph(builder);

  nc::net::ExclusionSet to_exclude;
  for (nc::net::GraphLinkIndex link : full_graph.AllLinks()) {
    nc::net::Delay delay = full_graph.GetLink(link)->delay();
    if (delay > threshold) {
      to_exclude.Links({link});
    }
  }

  std::vector<nc::net::GraphNodeSet> connected_components;
  nc::net::GraphNodeSet all_nodes = full_graph.AllNodes();
  for (nc::net::GraphNodeIndex node : all_nodes) {
    if (IsAlreadySeen(connected_components, node)) {
      continue;
    }

    nc::net::GraphNodeSet reachable_from_node =
        nc::net::ReachableNodes(node, full_graph, to_exclude);
    connected_components.emplace_back(reachable_from_node);
  }

  std::vector<nc::net::GraphBuilder> out;
  for (const nc::net::GraphNodeSet& connected_component :
       connected_components) {
    nc::net::GraphLinkSet links_in_component =
        LinksBetweenNodes(full_graph, connected_component);
    if (links_in_component.Empty()) {
      continue;
    }

    nc::net::GraphBuilder component_builder;
    for (nc::net::GraphLinkIndex link : links_in_component) {
      const nc::net::GraphLink* link_ptr = full_graph.GetLink(link);
      component_builder.AddLink({link_ptr->src_id(), link_ptr->dst_id(),
                                 link_ptr->bandwidth(), link_ptr->delay()});
    }
    out.emplace_back(component_builder);
  }

  return out;
}

static double GetDatapointForTopology(const nc::net::GraphStorage& graph) {
  return ctr::GetFractionOfPairsAboveLinkFraction(graph, FLAGS_sp_fraction,
                                                  FLAGS_link_fraction_limit);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<ctr::TopologyAndFilename> topologies = ctr::GetTopologyInputs();
  if (FLAGS_dump_stats) {
    DumpStats(topologies);
    return 0;
  }

  for (const ctr::TopologyAndFilename& topology_and_filename : topologies) {
    const std::string& topology_file = topology_and_filename.file;
    const nc::net::GraphStorage& graph = *topology_and_filename.graph;
    nc::net::GraphBuilder builder = graph.ToBuilder();

    if (FLAGS_dump_llpd) {
      LOG(INFO) << "Topology " << topology_file << " LLPD "
                << GetDatapointForTopology(graph) << " CD "
                << ctr::CapacityDiversity(graph);
      continue;
    }

    if (FLAGS_split_threshold_ms > 0) {
      nc::net::Delay threshold = std::chrono::duration_cast<nc::net::Delay>(
          std::chrono::milliseconds(FLAGS_split_threshold_ms));
      std::vector<nc::net::GraphBuilder> connected_components =
          SplitTopology(builder, threshold);
      CHECK(!connected_components.empty());
      if (connected_components.size() == 1) {
        continue;
      }

      for (uint32_t i = 0; i < connected_components.size(); ++i) {
        std::string serialized = connected_components[i].ToRepetita();
        std::string top_name = nc::StrCat(topology_file, "_cc_", i);
        nc::File::WriteStringToFileOrDie(serialized, top_name);
        LOG(INFO) << "Wrote CC of " << topology_file << " to " << top_name;
      }
    } else {
      builder.ScaleCapacity(FLAGS_link_capacity_scale);
      builder.ScaleDelay(FLAGS_link_delay_scale);

      std::string serialized = builder.ToRepetita(graph.NodeOrderOrDie());
      nc::File::WriteStringToFileOrDie(serialized, topology_file);
      LOG(INFO) << "Overwrote " << topology_file << " capacity scale "
                << FLAGS_link_capacity_scale << " delay scale "
                << FLAGS_link_delay_scale;
    }
  }
}
