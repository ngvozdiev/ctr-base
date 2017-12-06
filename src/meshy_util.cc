#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <functional>
#include <map>
#include <memory>
#include <ratio>
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
#include "ncode_common/src/viz/graph.h"
#include "ncode_common/src/viz/grapher.h"
#include "ncode_common/src/viz/web_page.h"
#include "common.h"
#include "opt/oversubscription_model.h"

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

static void DumpGraphToHTML(const std::string& out,
                            const nc::net::GraphStorage& graph) {
  std::vector<nc::viz::EdgeData> edges;
  for (nc::net::GraphLinkIndex link : graph.AllLinks()) {
    const nc::net::GraphLink* link_ptr = graph.GetLink(link);
    double delay_ms =
        std::chrono::duration<double, std::milli>(link_ptr->delay()).count();

    std::vector<double> loads = {1.0};
    edges.emplace_back(link, loads, link_ptr->ToStringNoPorts(), delay_ms);
  }

  nc::viz::HtmlPage page;
  nc::viz::DisplayMode display_mode("default");
  nc::viz::GraphToHTML(edges, {}, {display_mode}, graph, &page);
  nc::File::WriteStringToFile(page.Construct(), out);
}

static double TreenessFraction(const nc::net::GraphStorage& graph) {
  // Will assume all links are bidirectional.
  double links = graph.LinkCount();
  links = links / 2;

  double links_to_remove = (links - (graph.NodeCount() - 1));
  return links_to_remove / links;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<std::string> topology_files = GetTopologyFiles();
  CHECK(!topology_files.empty());

  std::vector<std::pair<double, std::string>> treeness_and_topologies;
  std::vector<double> treeness_fractions;
  for (const std::string& topology_file : topology_files) {
    nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(topology_file));
    builder.RemoveMultipleLinks();
    nc::net::GraphStorage graph(builder);

    if (graph.NodeCount() < 10) {
      LOG(INFO) << "Skipping " << topology_file;
      continue;
    }

    LOG(INFO) << "Processing " << topology_file << " size "
              << graph.NodeCount();
    double tf = TreenessFraction(graph);
    treeness_fractions.emplace_back(tf);
    treeness_and_topologies.emplace_back(tf, topology_file);

    DumpGraphToHTML(
        nc::StrCat("meshy_out/", nc::File::ExtractFileName(topology_file), "_",
                   tf, ".html"),
        graph);
  }

  std::sort(treeness_and_topologies.begin(), treeness_and_topologies.end());
  for (const auto& treeness_and_graph : treeness_and_topologies) {
    LOG(INFO) << treeness_and_graph.first << " at "
              << treeness_and_graph.second;
  }

  nc::viz::CDFPlot plot;
  plot.AddData("treeness", treeness_fractions);
  plot.PlotToDir("meshy_out");
}
