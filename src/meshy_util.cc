#include <gflags/gflags.h>
#include <stddef.h>
#include <chrono>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/mc_flow.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/algorithm.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/perfect_hash.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/viz/grapher.h"
#include "common.h"
#include "opt/oversubscription_model.h"

DEFINE_string(topology_files, "", "Topology files");

using NodePair = std::pair<nc::net::GraphNodeIndex, nc::net::GraphNodeIndex>;

static constexpr nc::net::Bandwidth kDefaultDemand =
    nc::net::Bandwidth::FromGBitsPerSecond(100);
static constexpr size_t kDefaultFlowCount = 1000;
static constexpr char kSuperSourceId[] = "SuperSource";
static constexpr char kSuperSinkId[] = "SuperSink";

static std::vector<std::string> GetTopologyFiles() {
  std::vector<std::string> out;
  std::vector<std::string> split = nc::Split(FLAGS_topology_files, ",");
  for (const std::string& piece : split) {
    std::vector<std::string> files = nc::Glob(piece);
    out.insert(out.end(), files.begin(), files.end());
  }

  return out;
}

// Routes a set of aggregates on their shortest paths and returns the set of
// links that bottleneck first when all aggregates are grown at the same rate.
static std::pair<nc::net::GraphLinkSet, nc::net::Bandwidth>
BottleneckLinksForAggregates(
    const std::set<NodePair>& aggregates, const nc::net::GraphStorage& graph,
    const std::map<NodePair, std::unique_ptr<nc::net::Walk>>& shortest_paths) {
  // Will create a fake traffic matrix and routing configuration that only
  // contains the shortest path for each aggregate. Will then use the
  // oversubscription model to figure out what links will be bottlenecked.
  ctr::TrafficMatrix tm(&graph);
  for (const auto& src_and_dst : aggregates) {
    tm.AddDemand(ctr::AggregateId(src_and_dst),
                 ctr::DemandAndFlowCount(kDefaultDemand, kDefaultFlowCount));
  }

  ctr::RoutingConfiguration routing_config(tm);
  for (const auto& src_and_dst : aggregates) {
    const nc::net::Walk& sp =
        *(nc::FindOrDieNoPrint(shortest_paths, src_and_dst));
    routing_config.AddRouteAndFraction(ctr::AggregateId(src_and_dst),
                                       {ctr::RouteAndFraction(&sp, 1.0)});
  }

  ctr::OverSubModel oversub_model(routing_config);
  const nc::net::GraphLinkMap<double>& links_to_load =
      oversub_model.link_to_load();

  nc::net::GraphLinkSet out;
  for (const auto& link_and_load : links_to_load) {
    double load = *(link_and_load.second);
    if (load > 0.99) {
      out.Insert(link_and_load.first);
    }
  }

  // Need to also figure out what the total flow through the network is.
  nc::net::Bandwidth total;
  const std::map<const nc::net::Walk*, nc::net::Bandwidth>& per_flow_bw_map =
      oversub_model.per_flow_bandwidth_map();
  for (const auto& path_and_bw : per_flow_bw_map) {
    nc::net::Bandwidth bandwidth = path_and_bw.second * kDefaultFlowCount;
    total += bandwidth;
  }

  return {out, total};
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

// Returns the flow between the super-source that connects the sources of a
// number of aggregates and another sources and the super-sink that connects the
// destinations of a number of aggregates and another destination.
static nc::net::Bandwidth PairwiseFlow(const std::set<NodePair>& aggregates,
                                       nc::net::GraphNodeIndex src,
                                       nc::net::GraphNodeIndex dst,
                                       const nc::net::GraphStorage& graph) {
  // Will first create super source / sinks.
  nc::net::GraphBuilder graph_builder = graph.ToBuilder();

  for (const auto& src_and_dst : aggregates) {
    const std::string& src_id = graph.GetNode(src_and_dst.first)->id();
    const std::string& dst_id = graph.GetNode(src_and_dst.second)->id();
    graph_builder.AddLink(
        {kSuperSourceId, src_id, nc::net::Bandwidth::Max(), nc::net::Delay(1)});
    graph_builder.AddLink(
        {dst_id, kSuperSinkId, nc::net::Bandwidth::Max(), nc::net::Delay(1)});
  }

  const std::string& src_id = graph.GetNode(src)->id();
  const std::string& dst_id = graph.GetNode(dst)->id();
  graph_builder.AddLink(
      {kSuperSourceId, src_id, nc::net::Bandwidth::Max(), nc::net::Delay(1)});
  graph_builder.AddLink(
      {dst_id, kSuperSinkId, nc::net::Bandwidth::Max(), nc::net::Delay(1)});

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
    const std::set<NodePair>& aggregates, const nc::net::GraphStorage& graph,
    const std::map<NodePair, std::unique_ptr<nc::net::Walk>>& shortest_paths) {
  double total = 0;
  double unreachable = 0;
  double unreachable_no_alternative = 0;

  nc::net::GraphLinkSet full_links;
  nc::net::Bandwidth total_flow;
  std::tie(full_links, total_flow) =
      BottleneckLinksForAggregates(aggregates, graph, shortest_paths);

  nc::net::ExclusionSet exclusion_set;
  exclusion_set.Links(full_links);

  nc::net::GraphNodeSet all_nodes = graph.AllNodes();
  for (nc::net::GraphNodeIndex src : all_nodes) {
    for (nc::net::GraphNodeIndex dst : all_nodes) {
      if (src == dst) {
        continue;
      }

      if (nc::ContainsKey(aggregates, std::make_pair(src, dst))) {
        continue;
      }

      ++total;
      if (Reachable(src, dst, graph, exclusion_set)) {
        continue;
      }

      ++unreachable;
      nc::net::Bandwidth pairwise_flow =
          PairwiseFlow(aggregates, src, dst, graph);
      if (pairwise_flow == total_flow) {
        ++unreachable_no_alternative;
      }
    }
  }

  return {unreachable / total, unreachable_no_alternative / total};
}

// Returns the aggregates whose shortest paths cross a link for each link.
nc::net::GraphLinkMap<std::set<NodePair>> SPCrossAggregates(
    const nc::net::GraphStorage& graph,
    std::map<NodePair, std::unique_ptr<nc::net::Walk>>* shortest_paths) {
  nc::net::GraphLinkMap<std::set<NodePair>> out;

  for (nc::net::GraphNodeIndex src : graph.AllNodes()) {
    nc::net::ShortestPath sp_tree(
        src, graph.AllNodes(), nc::net::ExclusionSet(), graph.AdjacencyList());
    for (nc::net::GraphNodeIndex dst : graph.AllNodes()) {
      if (src == dst) {
        continue;
      }

      std::unique_ptr<nc::net::Walk> shortest_path = sp_tree.GetPath(dst);
      for (nc::net::GraphLinkIndex link_on_path : shortest_path->links()) {
        out[link_on_path].emplace(src, dst);
      }

      (*shortest_paths)[{src, dst}] = std::move(shortest_path);
    }
  }

  return out;
}

static double PlotReachableFractions(const nc::net::GraphStorage& graph,
                                     const std::string& graph_name) {
  std::map<NodePair, std::unique_ptr<nc::net::Walk>> shortest_paths;
  nc::net::GraphLinkMap<std::set<NodePair>> cross_map =
      SPCrossAggregates(graph, &shortest_paths);

  std::vector<double> unreachable_values;
  std::vector<double> unreachable_no_alt_values;
  double total_delta = 0;

  for (const auto& link_and_aggregates : cross_map) {
    const std::set<NodePair> aggregates = *(link_and_aggregates.second);

    double unreachable;
    double unreachable_no_alternative;
    std::tie(unreachable, unreachable_no_alternative) =
        ReachableFraction(aggregates, graph, shortest_paths);

    unreachable_values.emplace_back(unreachable);
    unreachable_no_alt_values.emplace_back(unreachable_no_alternative);

    total_delta += unreachable - unreachable_no_alternative;
  }

  nc::viz::PythonGrapher grapher(nc::StrCat("meshy_out/g_", graph_name));
  grapher.PlotCDF({}, {{"unreachable", unreachable_values},
                       {"unreachable_no_alt", unreachable_no_alt_values}});
  double weight = graph.NodeCount() * (graph.NodeCount() - 1) / 2.0;
  return total_delta / weight;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<std::string> topology_files = GetTopologyFiles();
  CHECK(!topology_files.empty());

  std::vector<std::pair<double, std::string>> deltas_and_topologies;
  for (const std::string& topology_file : topology_files) {
    nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(topology_file));
    builder.RemoveMultipleLinks();
    nc::net::GraphStorage graph(builder);

    if (graph.NodeCount() > 70 || graph.NodeCount() < 20) {
      LOG(INFO) << "Skipping " << topology_file;
      continue;
    }

    LOG(INFO) << "Processing " << topology_file << " size "
              << graph.NodeCount();
    std::string filename = nc::File::ExtractFileName(topology_file);
    double delta = PlotReachableFractions(graph, filename);
    deltas_and_topologies.emplace_back(delta, filename);
  }

  std::sort(deltas_and_topologies.begin(), deltas_and_topologies.end());
  LOG(INFO) << nc::Join(
      deltas_and_topologies, "\n",
      [](const std::pair<double, std::string>& delta_and_topology) {
        return nc::StrCat(delta_and_topology.first, " : ",
                          delta_and_topology.second);
      });
}
