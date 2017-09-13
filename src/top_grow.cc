// Generates a series of topologies from a set of locations. Each of
// the topologies will have a progressively higher level of connectivity,
// starting from no connectivity and ending at a full clique. Each link will
// have a delay set up to be the speed of light in fiber.

#include <gflags/gflags.h>
#include <stddef.h>
#include <stdint.h>
#include <algorithm>
#include <chrono>
#include <limits>
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
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/thread_runner.h"
#include "ncode_common/src/viz/web_page.h"
#include "common.h"
#include "geo/geo.h"
#include "metrics/metrics.h"
#include "opt/ctr.h"
#include "opt/opt.h"
#include "opt/path_provider.h"

DEFINE_string(cities_file, "cities5000.txt", "Location of the cities file");
DEFINE_string(topology_file, "", "Topology file");
DEFINE_string(tm_file, "", "Traffic matrix file");
DEFINE_double(link_capacity_scale, 1.0, "By how much to scale all links");
DEFINE_double(delay_scale, 1.0, "By how much to scale the delays of all links");
DEFINE_double(
    light_speed, 200.0,
    "Speed of light (in km per millisecond), default is speed in fiber");
DEFINE_double(link_speed_Mbps, 1000, "All new links will be this fast");
DEFINE_string(
    optimizer, "CTR",
    "The optimizer to use. One of SP,CTR,MinMax,MinMaxLD,MinMaxK10,B4.");

using namespace std::chrono;

// Elements of a traffic matrix.
using TMElements = std::vector<nc::lp::DemandMatrixElement>;

// Number of flows in the most loaded aggregate.
static constexpr size_t kTopAggregateFlowCount = 10000;

static auto* total_delay_ms =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint32_t>("total_delay_ms", "Sum of all delays");

static auto* max_link_utilization =
    nc::metrics::DefaultMetricManager() -> GetThreadSafeMetric<double>(
        "max_link_utilization", "Maximum link utilization");

static auto* capacity_added =
    nc::metrics::DefaultMetricManager() -> GetThreadSafeMetric<double>(
        "capacity_added", "Capacity added to the network (in Mbps)");

namespace ctr {

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

static std::unique_ptr<RoutingConfiguration> RunOptimizer(
    const nc::lp::DemandMatrix& demand_matrix, PathProvider* provider) {
  auto tm = FromDemandMatrix(demand_matrix);

  std::unique_ptr<Optimizer> opt;
  if (FLAGS_optimizer == "CTR") {
    opt = nc::make_unique<CTROptimizer>(provider, 1.0, false, false);
  } else if (FLAGS_optimizer == "MinMax") {
    opt = nc::make_unique<MinMaxOptimizer>(provider, 1.0, false);
  } else if (FLAGS_optimizer == "MinMaxLD") {
    opt = nc::make_unique<MinMaxOptimizer>(provider, 1.0, true);
  } else if (FLAGS_optimizer == "MinMaxK10") {
    opt = nc::make_unique<MinMaxPathBasedOptimizer>(provider, 1.0, true, 10);
  } else if (FLAGS_optimizer == "B4") {
    opt = nc::make_unique<B4Optimizer>(provider, false, 1.0);
  } else if (FLAGS_optimizer == "SP") {
    opt = nc::make_unique<ShortestPathOptimizer>(provider);
  } else {
    LOG(FATAL) << "Bad optimizer " << FLAGS_optimizer;
  }

  return opt->Optimize(*tm);
}

// Scores a combination of topology and TM.
static std::pair<double, bool> Cost(const nc::net::GraphBuilder& topology,
                                    const TMElements& tm_elements) {
  nc::net::GraphStorage graph(topology);
  PathProvider path_provider(&graph);
  nc::lp::DemandMatrix demand_matrix(tm_elements, &graph);
  auto routing = RunOptimizer(demand_matrix, &path_provider);

  double total_oversubscription = 0;
  double max_utilization = 0;
  nc::net::GraphLinkMap<double> link_utilizations = routing->LinkUtilizations();
  for (const auto& link_and_utilization : link_utilizations) {
    double utilization = *(link_and_utilization.second);

    double oversubscription = utilization - 1.0;
    oversubscription = std::max(0.0, oversubscription);
    total_oversubscription += oversubscription;

    max_utilization = std::max(max_utilization, utilization);
  }

  return {total_oversubscription, max_utilization <= 1.0};
}

// Records information about the routing.
static void RecordRoutingConfig(const RoutingConfiguration& routing,
                                uint32_t step) {
  auto* total_path_handle = total_delay_ms->GetHandle();
  auto* link_utilization_handle = max_link_utilization->GetHandle();

  nc::net::Delay total_delay = routing.TotalPerFlowDelay();
  total_path_handle->AddValue(duration_cast<milliseconds>(total_delay).count());

  double max_utilization = routing.MaxLinkUtilization();
  link_utilization_handle->AddValue(max_utilization);

  nc::viz::HtmlPage page;
  routing.ToHTML(&page);
  nc::File::WriteStringToFile(page.Construct(),
                              nc::StrCat("top_grow_step_", step, ".html"));
}

// Will add a single link (or increase the capacity of existing one) that
// minimizes the overall cost.
static nc::net::GraphBuilder GrowNetwork(
    const nc::net::GraphBuilder& topology, const TMElements& tm_elements,
    const std::map<std::pair<std::string, std::string>, milliseconds>& delays) {
  nc::net::GraphBuilder to_return = topology;
  double total_added = 0.0;
  while (true) {
    nc::net::GraphBuilder best_topology;
    double best_cost = std::numeric_limits<double>::max();
    bool best_fits = false;
    std::string best_link_to_string;

    std::vector<nc::net::GraphBuilder> to_run;
    std::vector<std::string> links_added_to_string;
    std::set<std::string> all_nodes = topology.AllNodeNames();
    for (const std::string& src : all_nodes) {
      for (const std::string& dst : all_nodes) {
        if (src >= dst) {
          continue;
        }

        nc::net::GraphBuilder new_topology = to_return;
        nc::net::Delay delay =
            nc::FindOrDieNoPrint(delays, std::make_pair(src, dst));
        nc::net::Bandwidth bw =
            nc::net::Bandwidth::FromMBitsPerSecond(FLAGS_link_speed_Mbps);
        new_topology.AddLink({src, dst, bw, delay});
        new_topology.AddLink({dst, src, bw, delay});
        new_topology.RemoveMultipleLinks();
        to_run.emplace_back(new_topology);
        links_added_to_string.emplace_back(nc::StrCat(src, "->", dst));
      }
    }

    std::vector<std::unique_ptr<std::pair<double, bool>>> output =
        nc::RunInParallelWithResult<nc::net::GraphBuilder,
                                    std::pair<double, bool>>(
            to_run, [&tm_elements](const nc::net::GraphBuilder& builder) {
              return nc::make_unique<std::pair<double, bool>>(
                  Cost(builder, tm_elements));
            }, FLAGS_threads);

    for (size_t i = 0; i < to_run.size(); ++i) {
      const nc::net::GraphBuilder& new_topology = to_run[i];
      double cost;
      bool fits;
      std::tie(cost, fits) = *(output[i]);

      if (cost < best_cost) {
        best_topology = new_topology;
        best_cost = cost;
        best_link_to_string = links_added_to_string[i];
        best_fits = fits;
      }
    }

    CHECK(best_cost != std::numeric_limits<double>::max());
    LOG(INFO) << "Growing network, best cost " << best_cost << " will add "
              << best_link_to_string << " fits " << best_fits;
    total_added += FLAGS_link_speed_Mbps;

    to_return = best_topology;
    if (best_fits) {
      break;
    }
  }

  capacity_added->GetHandle()->AddValue(total_added);
  return to_return;
}

static void Run(
    const nc::net::GraphBuilder& topology, const TMElements& tm_elements,
    const std::map<std::pair<std::string, std::string>, milliseconds>& delays) {
  nc::net::GraphBuilder current_topology = topology;
  TMElements current_elements = tm_elements;

  for (size_t i = 0; i < FLAGS_steps; ++i) {
    // Will first run to record information.
    nc::net::GraphStorage graph(current_topology);
    PathProvider path_provider(&graph);

    nc::lp::DemandMatrix demand_matrix(current_elements, &graph);
    auto output = RunOptimizer(demand_matrix, &path_provider);
    RecordRoutingConfig(*output, i);
    LOG(INFO) << "Step " << i << " " << output->TotalPerFlowDelay().count()
              << " links " << current_topology.links().size();

    // Will grow the TM to figure out if any new links need to be added.
    auto scaled_demands = demand_matrix.Scale(FLAGS_capacity_grow_threshold);
    output = RunOptimizer(*scaled_demands, &path_provider);

    bool all_fit = output->MaxLinkUtilization() <= 1.0;
    if (!all_fit) {
      // Need to grow the network.
      current_topology =
          GrowNetwork(current_topology, scaled_demands->elements(), delays);
    } else {
      capacity_added->GetHandle()->AddValue(0.0);
    }

    // Perform a single step.
    scaled_demands = demand_matrix.Scale(FLAGS_scale_step);
    current_elements = scaled_demands->elements();
  }
}

}  // namespace ctr

static const nc::geo::CityData* GetCity(const std::string& name,
                                        nc::geo::Localizer* localizer) {
  nc::geo::FindCityRequest request;
  request.ascii_name = name;

  const nc::geo::CityData* city_data = localizer->FindCityOrNull(request);
  LOG(INFO) << name << " localized to " << city_data->ToString();
  CHECK(city_data != nullptr);
  return city_data;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::metrics::InitMetrics();
  nc::geo::Localizer localizer(FLAGS_cities_file);

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology_file), &node_order);
  builder.RemoveMultipleLinks();
  builder.ScaleCapacity(FLAGS_link_capacity_scale);
  builder.ScaleDelay(FLAGS_delay_scale);

  nc::net::GraphStorage graph(builder);
  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
      nc::lp::DemandMatrix::LoadRepetitaOrDie(
          nc::File::ReadFileToStringOrDie(FLAGS_tm_file), node_order, &graph);

  // Will record the delays between all N * (N - 1) nodes.
  std::map<std::pair<std::string, std::string>, milliseconds> delays;
  nc::net::GraphNodeSet all_nodes = graph.AllNodes();
  for (nc::net::GraphNodeIndex src : all_nodes) {
    for (nc::net::GraphNodeIndex dst : all_nodes) {
      if (src == dst) {
        continue;
      }

      const std::string& src_id = graph.GetNode(src)->id();
      const std::string& dst_id = graph.GetNode(dst)->id();

      const nc::geo::CityData* src_city = GetCity(src_id, &localizer);
      const nc::geo::CityData* dst_city = GetCity(dst_id, &localizer);
      double distance_km = src_city->DistanceKm(*dst_city);

      uint32_t delay_ms = distance_km / FLAGS_light_speed;
      delay_ms = std::max(delay_ms, static_cast<uint32_t>(1));
      delays[{src_id, dst_id}] = milliseconds(delay_ms);
    }
  }

  ctr::Run(builder, demand_matrix->elements(), delays);
}
