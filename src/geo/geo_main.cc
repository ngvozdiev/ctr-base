#include <gflags/gflags.h>
#include <string>

#include "ncode/file.h"
#include "ncode/logging.h"
#include "ncode/net/net_common.h"
#include "ncode/net/net_gen.h"
#include "geo.h"

DEFINE_string(cities_file, "cities5000.txt", "Location of the cities file");
DEFINE_string(topology_file, "", "Location of the topology to check");
DEFINE_string(output_file, "", "Where to store the new topology");

using namespace nc::geo;
using namespace std::chrono;

static const CityData* Lookup(const std::string& id, Localizer* localizer) {
  FindCityRequest src_request;
  src_request.ascii_name = id;
  return localizer->FindCityOrNull(src_request);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_topology_file.empty());

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology_file), &node_order);
  builder.RemoveMultipleLinks();
  nc::net::GraphStorage graph(builder);

  Localizer localizer(FLAGS_cities_file);
  WebMercator projection;

  nc::net::GraphBuilder new_builder;
  for (nc::net::GraphLinkIndex link : graph.AllLinks()) {
    const nc::net::GraphLink* link_ptr = graph.GetLink(link);
    const std::string& src = link_ptr->src_id();
    const std::string& dst = link_ptr->dst_id();

    const CityData* src_city = Lookup(src, &localizer);
    const CityData* dst_city = Lookup(dst, &localizer);
    if (src_city == nullptr || dst_city == nullptr) {
      LOG(INFO) << "Unable to find " << link_ptr->ToStringNoPorts();
      continue;
    }

    milliseconds link_delay = duration_cast<milliseconds>(link_ptr->delay());

    double distance = src_city->DistanceKm(*dst_city);

    // Light travels 200km per millisecond in fiber.
    double delay_ms = distance / 200.0;

    LOG(INFO) << "Link from " << src << " (" << src_city->name << ") to " << dst
              << " (" << dst_city->name << ")"
              << " distance " << distance << "km speed of light latency "
              << delay_ms << "ms vs " << link_delay.count() << "ms";

    uint64_t delay_int = static_cast<uint64_t>(delay_ms);
    delay_int = std::max(static_cast<uint64_t>(1), delay_int);

    new_builder.AddLink({src, dst, link_ptr->bandwidth(),
                         std::chrono::milliseconds(delay_int)});
  }

  if (!FLAGS_output_file.empty()) {
    nc::File::WriteStringToFileOrDie(new_builder.ToRepetita(node_order),
                                     FLAGS_output_file);
  }
}
