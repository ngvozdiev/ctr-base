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
DEFINE_double(delay_scale, 1.0, "By how much to scale the delays of all links");

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

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<std::string> topology_files = GetTopologyFiles();
  CHECK(!topology_files.empty());

  std::vector<double> fractions;
  for (const std::string& topology_file : topology_files) {
    nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(topology_file));
    builder.RemoveMultipleLinks();
    builder.ScaleDelay(FLAGS_delay_scale);
    nc::net::GraphStorage graph(builder);

    double common_fraction;
    std::tie(std::ignore, common_fraction) = nc::net::CommonSPLinks(graph);
    fractions.emplace_back(common_fraction);

    LOG(INFO) << topology_file << " " << common_fraction << " "
              << graph.AllNodes().Count();

    if (topology_files.size() == 1) {
      DumpGraphToHTML("out.html", graph);
    }
  }

  if (topology_files.size() > 1) {
    nc::viz::DataSeries1D data_1d;
    data_1d.data = std::move(fractions);

    nc::viz::PythonGrapher python_grapher("mesy_fractions");
    python_grapher.PlotCDF({}, {data_1d});
  }
}
