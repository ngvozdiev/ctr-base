#include <gflags/gflags.h>
#include <chrono>
#include <iostream>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/mc_flow.h"
#include "ncode_common/src/net/algorithm.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/perfect_hash.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/viz/grapher.h"

DEFINE_string(topology_files, "", "Topology files");

static std::vector<std::string> GetTopologyFiles() {
  std::vector<std::string> out;
  std::vector<std::string> split = nc::Split(FLAGS_topology_files, ",");
  for (const std::string& piece : split) {
    std::vector<std::string> files = nc::Glob(piece);
    out.insert(out.end(), files.begin(), files.end());
  }

  return out;
}

// Returns true if two nodes are reachable given a set of constraints.
static bool Reachable(nc::net::GraphNodeIndex node_one,
                      nc::net::GraphNodeIndex node_two,
                      const nc::net::GraphStorage& graph,
                      const nc::net::ExclusionSet& exclusion_set) {
  nc::net::GraphNodeSet reachable_from_one =
      nc::net::ReachableNodes(node_one, graph, exclusion_set);
  return reachable_from_one.Contains(node_two);
}

static constexpr char kSuperSourceId[] = "SuperSource";
static constexpr char kSuperSinkId[] = "SuperSink";

// Returns the flow between the super-source that connects two sources and the
// super-sink that connects two sink nodes.
static nc::net::Bandwidth PairwiseFlow(nc::net::GraphNodeIndex source_one,
                                       nc::net::GraphNodeIndex source_two,
                                       nc::net::GraphNodeIndex sink_one,
                                       nc::net::GraphNodeIndex sink_two,
                                       const nc::net::GraphStorage& graph) {
  const std::string& source_one_id = graph.GetNode(source_one)->id();
  const std::string& source_two_id = graph.GetNode(source_two)->id();
  const std::string& sink_one_id = graph.GetNode(sink_one)->id();
  const std::string& sink_two_id = graph.GetNode(sink_two)->id();

  // Will first create super source / sinks.
  nc::net::GraphBuilder graph_builder = graph.ToBuilder();
  graph_builder.AddLink({kSuperSourceId, source_one_id,
                         nc::net::Bandwidth::Max(), nc::net::Delay(1)});
  graph_builder.AddLink({kSuperSourceId, source_two_id,
                         nc::net::Bandwidth::Max(), nc::net::Delay(1)});
  graph_builder.AddLink({sink_one_id, kSuperSinkId, nc::net::Bandwidth::Max(),
                         nc::net::Delay(1)});
  graph_builder.AddLink({sink_two_id, kSuperSinkId, nc::net::Bandwidth::Max(),
                         nc::net::Delay(1)});

  // The flow problem will take link capacities in Mbps.
  nc::net::GraphStorage extended_graph(graph_builder);
  nc::net::GraphLinkMap<double> link_capacities;
  for (nc::net::GraphLinkIndex link : extended_graph.AllLinks()) {
    const nc::net::GraphLink* link_ptr = extended_graph.GetLink(link);
    double capacity = link_ptr->bandwidth().Mbps();
    link_capacities[link] = capacity;
  }
  nc::lp::MaxFlowSingleCommodityFlowProblem flow_problem(link_capacities,
                                                         &extended_graph);
  flow_problem.AddDemand(kSuperSourceId, kSuperSinkId);

  double max_flow_Mbps;
  CHECK(flow_problem.GetMaxFlow(&max_flow_Mbps));
  return nc::net::Bandwidth::FromMBitsPerSecond(max_flow_Mbps);
}

// Returns the fraction of all pairs that are reachable when the shortest path
// between two nodes in excluded from the graph.
static std::pair<double, double> ReachableFraction(
    nc::net::GraphNodeIndex node_one, nc::net::GraphNodeIndex node_two,
    const nc::net::GraphStorage& graph) {
  double total = 0;
  double unreachable = 0;
  double unreachable_no_alternative = 0;

  nc::net::ShortestPath sp_tree(node_one, graph.AllNodes(), {},
                                graph.AdjacencyList());
  std::unique_ptr<nc::net::Walk> path = sp_tree.GetPath(node_two);
  CHECK(path && path->size() != 0) << "Unable to find shortest path";

  nc::net::Bandwidth bottleneck_bandwidth;
  nc::net::ExclusionSet exclusion_set;
  exclusion_set.Links(path->BottleneckLinks(graph, &bottleneck_bandwidth));

  nc::net::GraphNodeSet all_nodes = graph.AllNodes();
  for (nc::net::GraphNodeIndex src : all_nodes) {
    for (nc::net::GraphNodeIndex dst : all_nodes) {
      if (src <= dst) {
        // Don't care about half of the possibilities, since links are
        // bidirectional.
        continue;
      }

      if (src == node_one && dst == node_two) {
        continue;
      }

      ++total;
      if (Reachable(src, dst, graph, exclusion_set)) {
        continue;
      }

      ++unreachable;
      nc::net::Bandwidth pairwise_flow =
          PairwiseFlow(node_one, src, node_two, dst, graph);
      if (pairwise_flow == bottleneck_bandwidth) {
        ++unreachable_no_alternative;
      }
    }
  }

  return {unreachable / total, unreachable_no_alternative / total};
}

static double PlotReachableFractions(const nc::net::GraphStorage& graph,
                                     const std::string& graph_name) {
  std::vector<double> unreachable_values;
  std::vector<double> unreachable_no_alt_values;
  double total_delta = 0;

  nc::net::GraphNodeSet all_nodes = graph.AllNodes();
  for (nc::net::GraphNodeIndex src : all_nodes) {
    for (nc::net::GraphNodeIndex dst : all_nodes) {
      if (src <= dst) {
        continue;
      }

      double unreachable;
      double unreachable_no_alt;
      std::tie(unreachable, unreachable_no_alt) =
          ReachableFraction(src, dst, graph);
      unreachable_values.emplace_back(unreachable);
      unreachable_no_alt_values.emplace_back(unreachable_no_alt);

      total_delta += unreachable - unreachable_no_alt;
    }
  }

  nc::viz::PythonGrapher grapher(nc::StrCat("meshy_out/g_", graph_name));
  grapher.PlotCDF({}, {{"unreachable", unreachable_values},
                       {"unreachable_no_alt", unreachable_no_alt_values}});
  double weight = all_nodes.Count() * (all_nodes.Count() - 1) / 2.0;
  return total_delta / weight;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<std::string> topology_files = GetTopologyFiles();
  CHECK(!topology_files.empty());

  double max_delta = 0;
  std::string topology_with_max_delta;

  std::vector<double> all_links_to_remove;
  for (const std::string& topology_file : topology_files) {
    nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(topology_file));
    builder.RemoveMultipleLinks();
    nc::net::GraphStorage graph(builder);

    if (graph.NodeCount() > 50 || graph.NodeCount() < 20) {
      LOG(INFO) << "Skipping " << topology_file;
      continue;
    }

    LOG(INFO) << "Processing " << topology_file << " size "
              << graph.NodeCount();
    std::string filename = nc::File::ExtractFileName(topology_file);
    double delta = PlotReachableFractions(graph, filename);

    if (delta > max_delta) {
      max_delta = delta;
      topology_with_max_delta = filename;
    }
  }

  LOG(INFO) << "max delta " << max_delta << " topology "
            << topology_with_max_delta;
}
