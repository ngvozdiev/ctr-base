// Generates a series of topologies from a set of locations. Each of
// the topologies will have a progressively higher level of connectivity,
// starting from no connectivity and ending at a full clique. Each link will
// have a delay set up to be the speed of light in fiber.

#include <gflags/gflags.h>
#include <stdint.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/substitute.h"
#include "geo/geo.h"

DEFINE_string(cities_file, "cities5000.txt", "Location of the cities file");
DEFINE_string(locations, "", "Locations (comma-separated)");
DEFINE_uint64(seed, 1, "Seed to use when growing the topology");
DEFINE_string(output, "top_grow_$0.graph",
              "Output, $0 will be replaced with the number of links in the "
              "topology (will range from 1 to len(location) * (len(locations) "
              "- 1))");
DEFINE_double(
    light_speed, 200.0,
    "Speed of light (in km per millisecond), default is speed in fiber");
DEFINE_double(link_speed_Mbps, 1000, "Speed of all links will be this value");

using namespace std::chrono;
using Endpoint = std::pair<const nc::geo::CityData*, uint32_t>;

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<std::string> locations = nc::Split(FLAGS_locations, ",");

  nc::geo::Localizer localizer(FLAGS_cities_file);
  nc::geo::WebMercator projection;

  std::vector<Endpoint> cities;
  for (const std::string& location : locations) {
    nc::geo::FindCityRequest request;
    request.ascii_name = location;

    const nc::geo::CityData* city_data = localizer.FindCityOrNull(request);
    LOG(INFO) << location << " localized to " << city_data->ToString();
    CHECK(city_data != nullptr);
    uint32_t i = cities.size();
    cities.push_back({city_data, i});
  }

  // Will randomly shuffle the cities to get a list of starting nodes, will then
  // for each one get another randomly shuffled list to get end nodes.
  std::vector<Endpoint> srcs = cities;
  std::mt19937 rnd(FLAGS_seed);
  std::shuffle(srcs.begin(), srcs.end(), rnd);

  nc::net::GraphBuilder builder;
  for (Endpoint src : srcs) {
    std::vector<Endpoint> dsts = cities;
    std::shuffle(dsts.begin(), dsts.end(), rnd);
    for (Endpoint dst : dsts) {
      if (src == dst) {
        continue;
      }

      const nc::geo::CityData* src_city_data = src.first;
      const nc::geo::CityData* dst_city_data = dst.first;
      double distance_km = src_city_data->DistanceKm(*dst_city_data);
      distance_km = std::max(distance_km, 200.0);
      uint32_t delay_micros = (distance_km / FLAGS_light_speed) * 1000.0;

      std::string src_name = nc::StrCat(src_city_data->name, src.second);
      std::string dst_name = nc::StrCat(dst_city_data->name, dst.second);

      builder.AddLink(
          {src_name, dst_name,
           nc::net::Bandwidth::FromMBitsPerSecond(FLAGS_link_speed_Mbps),
           duration_cast<nc::net::Delay>(microseconds(delay_micros))});
      size_t link_count = builder.links().size();
      std::string serialized_topology = builder.ToRepetita();
      nc::File::WriteStringToFileOrDie(
          serialized_topology,
          nc::Substitute(FLAGS_output.c_str(), link_count));
    }
  }
}
