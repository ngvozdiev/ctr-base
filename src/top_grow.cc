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
DEFINE_string(topology_file, "", "Topology file");
DEFINE_string(tm_file, "", "Traffic matrix file");
DEFINE_uint64(seed, 1, "Seed to use when growing the topology");
DEFINE_string(output, "top_grow_$0.graph",
              "Output, $0 will be replaced with the number of links in the "
              "topology (will range from 1 to len(location) * (len(locations) "
              "- 1))");
DEFINE_double(
    light_speed, 200.0,
    "Speed of light (in km per millisecond), default is speed in fiber");
DEFINE_double(link_speed_Mbps, 1000, "All new links will be this fast");
DEFINE_double(scale_step, 1.01, "By how much to scale the TM each step");
DEFINE_double(capacity_grow_threshold, 0.3,
              "What the limit is for growing capacity");

using namespace std::chrono;
using Endpoint = std::pair<const nc::geo::CityData*, uint32_t>;

static auto* total_delay_ms =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint32_t, uint32_t>("total_delay_ms",
                                                   "Sum of all delays", "Step");

namespace ctr {

// Elements of a traffic matrix.
using TMElements = std::vector<nc::lp::DemandMatrixElement>;

// Number of flows in the most loaded aggregate.
static constexpr size_t kTopAggregateFlowCount = 10000;

static std::unique_ptr<TrafficMatrix> FromDemandMatrix(
    const nc::lp::DemandMatrix& demand_matrix) {
  nc::net::Bandwidth max_demand = nc::net::Bandwidth::Zero();
  for (const auto& element : demand_matrix.elements()) {
    max_demand = std::max(max_demand, element.demand);
  }

  // Will assign flow counts so that the max bandwidth aggregate has
  // kTopAggregateFlowCount, and all other aggregates proportionally less.
  std::map<AggregateId, DemandAndFlowCount> demands_and_counts;
  for (const auto& element : demand_matrix.elements()) {
    size_t flow_count = kTopAggregateFlowCount * (element.demand / max_demand);
    flow_count = std::max(1ul, flow_count);

    AggregateId id(element.src, element.dst);
    demands_and_counts[id] = {element.demand, flow_count};
  }

  return nc::make_unique<TrafficMatrix>(demand_matrix.graph(),
                                        demands_and_counts);
}

std::unique_ptr<RoutingConfiguration> RunCTR(
    const nc::net::GraphStorage& graph,
    const nc::lp::DemandMatrix demand_matrix) {
  auto tm = FromDemandMatrix(demand_matrix);
  PathProvider path_provider(graph);
  CTROptimizer ctr_optimizer(&path_provider, 1.0, false, false);
  return ctr_optimizer.Optimize(*tm);
}

// Scores a combination of topology and TM.
static double Cost(const nc::net::GraphBuilder& topology,
                   const TMElements& tm_elements);

// Records information about the routing.
static void RecordRoutingConfig(const RoutingConfiguration& routing,
                                size_t step) {
  auto* total_path_handle = total_delay_ms->GetHandle(step);
  nc::net::Delay total_delay = routing.TotalPerFlowDelay();
  total_path_handle->AddValue(duration_cast<milliseconds>(total_delay).count());
}

// Runs the optimizer and records stuff.
static void RunOptimizer(const nc::net::GraphBuilder& topology,
                         const TMElements& tm_elements, size_t step) {
  nc::net::GraphStorage graph(topology);
  nc::lp::DemandMatrix demand_matrix(tm_elements, &graph);
  std::unique_ptr<RoutingConfiguration> routing = RunCTR(graph, demand_matrix);
  RecordRoutingConfig(*routing, step);
}

// Will add a single link (or increase the capacity of existing one) that
// minimizes the overall cost.
static nc::net::GraphBuilder GrowNetwork(
    const nc::net::GraphBuilder& topology, const TMElements& tm_elements,
    const std::map<std::pair<std::string, std::string>, milliseconds>& delays) {
  double best_cost = std::numeric_limits<double>::max();
  nc::net::GraphBuilder best_topology;

  std::set<std::string> all_nodes = topology.AllNodeNames();
  for (const std::string& src : all_nodes) {
    for (const std::string& dst : all_nodes) {
      if (src == dst) {
        continue;
      }

      nc::net::GraphBuilder new_topology = topology;
      nc::net::Delay delay = nc::FindOrDie(delays, {src, dst});
      nc::net::Bandwidth bw =
          nc::net::Bandwidth::FromMBitsPerSecond(FLAGS_link_speed_Mbps);
      new_topology.AddLink({src, dst, bw, delay});
      new_topology.RemoveMultipleLinks();

      double cost = Cost(new_topology, tm_elements);
      if (cost < best_cost) {
        best_topology = new_topology;
        best_cost = cost;
      }
    }
  }

  CHECK(best_cost != std::numeric_limits<double>::max());
  return best_topology;
}

static void Run(
    const nc::net::GraphBuilder& topology, const TMElements& tm_elements,
    const std::map<std::pair<std::string, std::string>, milliseconds>& delays) {
  nc::net::GraphBuilder current_topology = topology;
  TMElements current_elements = tm_elements;

  for (size_t i = 0; i < 100; ++i) {
    // Will first grow the TM.
    nc::net::GraphStorage graph(current_topology);
    nc::lp::DemandMatrix demand_matrix(current_elements, &graph);
    auto scaled_demands = demand_matrix.Scale(FLAGS_scale_step);

    nc::net::Bandwidth flow;
    double max_commodity_scale;
    std::tie(flow, max_commodity_scale) = scaled_demands->GetMaxFlow({});
    CHECK(flow != nc::net::Bandwidth::Zero());

    const nc::net::GraphBuilder* topology_ptr = &topology;
    std::unique_ptr<nc::net::GraphBuilder> new_topology;
    if (max_commodity_scale < FLAGS_capacity_grow_threshold) {
      // Time to grow the network.
      new_topology = nc::make_unique<nc::net::GraphBuilder>(
          GrowNetwork(topology, scaled_demands->elements(), delays));
      topology_ptr = new_topology.get();
    }

    // Run the optimizer with the new demands.
    RunOptimizer(*topology_ptr, scaled_demands->elements(), i);

    current_topology = *topology_ptr;
    current_elements = scaled_demands->elements();
  }
}

}  // namespace ctr

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
