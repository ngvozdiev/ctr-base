#include <gflags/gflags.h>
#include <ncode/common.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/map_util.h>
#include <ncode/net/net_common.h>
#include <ncode/net/net_gen.h>
#include <ncode/thread_runner.h>
#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "geo/geo.h"
#include "opt/opt.h"

DEFINE_string(topology_root, "", "Root for topologies. Required.");
DEFINE_double(fraction_links_to_add, 0.1, "Fraction of links to add");
DEFINE_double(sp_fraction, 1.4, "How far from the SP a path can be");
DEFINE_double(link_fraction_limit, 0.7,
              "At least this much of the SP's links can be routed around");

using namespace nc::geo;
using namespace std::chrono;

static constexpr char kTopologyExtension[] = ".graph";
static constexpr double kSpeedOfLightFiber = 392728119.98;

static std::vector<std::string> GetTopologyFiles() {
  std::string top_root = FLAGS_topology_root;
  CHECK(top_root != "");
  return nc::File::GetFilesWithExtension(top_root, kTopologyExtension);
}

static nc::net::Bandwidth MedianLinkCapacity(
    const nc::net::GraphStorage& graph) {
  std::vector<nc::net::Bandwidth> bw_values;
  for (nc::net::GraphLinkIndex link : graph.AllLinks()) {
    nc::net::Bandwidth bw = graph.GetLink(link)->bandwidth();
    bw_values.emplace_back(bw);
  }

  std::sort(bw_values.begin(), bw_values.end());
  CHECK(!bw_values.empty());
  return bw_values[bw_values.size() / 2];
}

const nc::geo::CityData* GetCityData(const std::string& name,
                                     nc::geo::Localizer* localizer) {
  nc::geo::FindCityRequest request;
  request.ascii_name = name;
  const nc::geo::CityData* city_data = localizer->FindCityOrNull(request);
  CHECK(city_data != nullptr);

  return city_data;
}

static std::unique_ptr<nc::net::GraphStorage> ExtendGraph(
    const nc::net::GraphStorage& graph, const std::string& src,
    const std::string& dst, nc::net::Bandwidth link_capacity,
    const std::map<std::string, const CityData*>& city_data) {
  const nc::geo::CityData* src_city = nc::FindOrDie(city_data, src);
  const nc::geo::CityData* dst_city = nc::FindOrDie(city_data, dst);
  double distance_m = src_city->DistanceKm(*dst_city) * 1000.0;
  double delay_s = distance_m / kSpeedOfLightFiber;
  double delay_ms = std::max(1.0, delay_s / 1000.0);
  nc::net::Delay delay = duration_cast<nc::net::Delay>(
      milliseconds(static_cast<size_t>(delay_ms)));

  nc::net::GraphBuilder builder = graph.ToBuilder();
  builder.AddLink({src, dst, link_capacity, delay});
  builder.AddLink({dst, src, link_capacity, delay});
  return nc::make_unique<nc::net::GraphStorage>(builder);
}

static std::unique_ptr<nc::net::GraphStorage> GrowTopology(
    const nc::net::GraphStorage& original_graph,
    nc::net::Bandwidth link_capacity,
    const std::map<std::string, const CityData*>& city_data) {
  LOG(INFO) << "Current " << ctr::GetFractionOfPairsAboveLinkFraction(
                                 original_graph, FLAGS_sp_fraction,
                                 FLAGS_link_fraction_limit);

  std::vector<std::pair<std::string, std::string>> inputs;
  for (nc::net::GraphNodeIndex src : original_graph.AllNodes()) {
    for (nc::net::GraphNodeIndex dst : original_graph.AllNodes()) {
      if (src == dst) {
        continue;
      }

      std::string src_id = original_graph.GetNode(src)->id();
      std::string dst_id = original_graph.GetNode(dst)->id();
      if (original_graph.HasLink(src_id, dst_id)) {
        continue;
      }

      if (src_id.find("_None_") != std::string::npos ||
          dst_id.find("_None_") != std::string::npos) {
        continue;
      }

      if (src_id.find("_?_") != std::string::npos ||
          dst_id.find("_?_") != std::string::npos) {
        continue;
      }

      if (src_id.find("_???_") != std::string::npos ||
          dst_id.find("_???_") != std::string::npos) {
        continue;
      }

      inputs.emplace_back(src_id, dst_id);
    }
  }

  std::vector<std::unique_ptr<double>> outputs =
      nc::RunInParallelWithResult<std::pair<std::string, std::string>, double>(
          inputs, [&original_graph, &link_capacity, &city_data](
                      const std::pair<std::string, std::string>& src_and_dst) {
            std::string src_id;
            std::string dst_id;
            std::tie(src_id, dst_id) = src_and_dst;

            auto new_graph = ExtendGraph(original_graph, src_id, dst_id,
                                         link_capacity, city_data);
            double llpd = ctr::GetFractionOfPairsAboveLinkFraction(
                *new_graph, FLAGS_sp_fraction, FLAGS_link_fraction_limit);
            LOG(INFO) << "L " << src_id << " -> " << dst_id << " " << llpd;
            return nc::make_unique<double>(llpd);
          });

  double max_llpd = 0;
  std::unique_ptr<nc::net::GraphStorage> max_graph;
  std::string max_src_id;
  std::string max_dst_id;
  for (size_t i = 0; i < outputs.size(); ++i) {
    const std::pair<std::string, std::string>& src_and_dst = inputs[i];
    std::string src_id;
    std::string dst_id;
    std::tie(src_id, dst_id) = src_and_dst;
    double llpd = *(outputs[i]);

    if (llpd > max_llpd) {
      auto new_graph =
          ExtendGraph(original_graph, src_id, dst_id, link_capacity, city_data);

      max_graph = std::move(new_graph);
      max_src_id = src_id;
      max_dst_id = dst_id;
      max_llpd = llpd;
    }
  }

  LOG(INFO) << "MAX " << max_src_id << " -> " << max_dst_id << " " << max_llpd;
  CHECK(max_graph);
  return max_graph;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::geo::Localizer localizer("cities5000.txt");

  std::vector<std::string> topology_files = GetTopologyFiles();
  for (const std::string& topology_file : topology_files) {
    std::vector<std::string> node_order;
    nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(topology_file), &node_order);

    auto graph = nc::make_unique<nc::net::GraphStorage>(builder);
    std::map<std::string, const CityData*> city_data;
    for (nc::net::GraphNodeIndex node : graph->AllNodes()) {
      const std::string id = graph->GetNode(node)->id();
      city_data[id] = GetCityData(id, &localizer);
    }

    size_t links_to_add =
        std::ceil(graph->LinkCount() * FLAGS_fraction_links_to_add);
    for (size_t i = 0; i < links_to_add; ++i) {
      graph = GrowTopology(*graph, MedianLinkCapacity(*graph), city_data);
    }

    nc::net::GraphBuilder new_builder = graph->ToBuilder();
    std::string serialized = new_builder.ToRepetita(node_order);
    std::string name = nc::File::ExtractFileName(topology_file);
    nc::File::WriteStringToFileOrDie(serialized,
                                     nc::StrCat("llpd_grow/", name));
  }
}
