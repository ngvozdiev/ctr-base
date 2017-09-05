#include <gflags/gflags.h>
#include <algorithm>
#include <chrono>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/algorithm.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/viz/graph.h"
#include "ncode_common/src/viz/grapher.h"
#include "ncode_common/src/viz/web_page.h"

DEFINE_string(topology_files, "", "Topology files");
DEFINE_double(delay_scale, 1.3, "By how much to scale the delays of all links");

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

// Returns the link with max delay.
const nc::net::GraphLink* MaxDelayLink(const nc::net::GraphStorage& graph) {
  nc::net::Delay max_delay = nc::net::Delay::zero();
  const nc::net::GraphLink* out = nullptr;
  for (const auto& link_index : graph.AllLinks()) {
    const nc::net::GraphLink* link = graph.GetLink(link_index);
    if (link->delay() > max_delay) {
      out = link;
      max_delay = link->delay();
    }
  }
  CHECK(out != nullptr);
  return out;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<std::string> topology_files = GetTopologyFiles();
  CHECK(!topology_files.empty());

  std::vector<double> all_links_to_remove;
  for (const std::string& topology_file : topology_files) {
    nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(topology_file));
    builder.RemoveMultipleLinks();
    builder.ScaleDelay(FLAGS_delay_scale);
    nc::net::GraphStorage graph(builder);

    const nc::net::GraphLink* max_delay_link = MaxDelayLink(graph);
    if (graph.AllNodes().Count() < 10) {
      continue;
    }

    size_t links_in_tree = graph.AllNodes().Count() - 1;
    size_t bidirectional_links = graph.AllLinks().Count() / 2;
    double links_to_remove = bidirectional_links - links_in_tree;
    all_links_to_remove.emplace_back(links_to_remove / bidirectional_links);

    LOG(INFO) << topology_file << " " << graph.AllNodes().Count() << " "
              << graph.AllLinks().Count() << " max delay link "
              << max_delay_link->ToStringNoPorts() << " delay "
              << max_delay_link->delay().count();

    if (topology_files.size() == 1) {
      DumpGraphToHTML("out.html", graph);
    }
  }

  if (topology_files.size() > 1) {
    nc::viz::DataSeries1D data_1d;
    data_1d.data = std::move(all_links_to_remove);

    nc::viz::PythonGrapher python_grapher("meshy_ltr");
    python_grapher.PlotCDF({}, {data_1d});
  }
}
