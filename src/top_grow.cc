#include <gflags/gflags.h>
#include <stddef.h>
#include <stdint.h>
#include <algorithm>
#include <chrono>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/viz/web_page.h"
#include "common.h"
#include "geo/geo.h"
#include "metrics/metrics.h"
#include "opt/ctr.h"
#include "opt/opt.h"
#include "opt/path_provider.h"
#include "opt_eval.h"

DEFINE_string(cities_file, "cities5000.txt", "Location of the cities file");
DEFINE_double(
    light_speed, 200.0,
    "Speed of light (in km per millisecond), default is speed in fiber");
DEFINE_double(link_speed_Mbps, 1000, "All new links will be this fast");
DEFINE_string(optimizers, "SP,CTR,CTRNFC,MinMax,MinMaxLD,MinMaxK10,B4",
              "The optimizers to use, comma-separated.");

using namespace std::chrono;

static auto* decreasing_delay_aggreggate =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string,
                               bool>(
            "decreasing_delay_aggregate_fraction",
            "Fraction of aggregates whose delay decreases", "Topology",
            "Traffic matrix", "Optimizer", "New link");

static auto* increasing_delay_aggregate =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string,
                               bool>(
            "increasing_delay_aggregate_fraction",
            "Fraction of aggregates whose delay increases", "Topology",
            "Traffic matrix", "Optimizer", "New link");

static auto* increasing_overload_delta =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string,
                               bool>(
            "overload_delta_fraction",
            "Fraction of aggregates that encounter more congestion", "Topology",
            "Traffic matrix", "Optimizer", "New link");

static auto* total_delay_delta =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string,
                               bool>("total_delay_delta_fraction",
                                     "Change in total delay", "Topology",
                                     "Traffic matrix", "Optimizer", "New link");

static auto* absolute_path_stretch_micros =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<nc::DiscreteDistribution<int64_t>, std::string,
                               std::string, std::string, bool>(
            "absolute_path_stretch_micros",
            "Distribution of absolute per-flow deltas in microseconds",
            "Topology", "Traffic matrix", "Optimizer", "New link");

// Elements of a traffic matrix.
using TMElements = std::vector<nc::lp::DemandMatrixElement>;

// Number of flows in the most loaded aggregate.
static constexpr size_t kTopAggregateFlowCount = 10000;

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
    const std::string& opt_string, const nc::lp::DemandMatrix& demand_matrix,
    PathProvider* provider) {
  auto tm = FromDemandMatrix(demand_matrix);

  std::unique_ptr<Optimizer> opt;
  if (opt_string == "CTR") {
    opt = nc::make_unique<CTROptimizer>(provider, 1.0, false, false);
  } else if (opt_string == "CTRNFC") {
    opt = nc::make_unique<CTROptimizer>(provider, 1.0, false, true);
  } else if (opt_string == "MinMax") {
    opt = nc::make_unique<MinMaxOptimizer>(provider, 1.0, false);
  } else if (opt_string == "MinMaxLD") {
    opt = nc::make_unique<MinMaxOptimizer>(provider, 1.0, true);
  } else if (opt_string == "MinMaxK10") {
    opt = nc::make_unique<MinMaxPathBasedOptimizer>(provider, 1.0, true, 10);
  } else if (opt_string == "B4") {
    opt = nc::make_unique<B4Optimizer>(provider, false, 1.0);
  } else if (opt_string == "SP") {
    opt = nc::make_unique<ShortestPathOptimizer>(provider);
  } else {
    LOG(FATAL) << "Bad optimizer " << opt_string;
  }

  return opt->Optimize(*tm);
}

static void RecordStretch(
    const OptEvalInput& input, const RoutingConfiguration& routing,
    const std::vector<const RoutingConfiguration*> new_routings,
    const std::string& opt, bool new_link) {
  nc::DiscreteDistribution<int64_t> dist;

  for (const RoutingConfiguration* new_routing : new_routings) {
    RoutingConfigurationDelta delta = routing.GetDifference(*new_routing);
    for (const auto& aggregate_and_delta : delta.aggregates) {
      // Need to figure out how many flows there are in an aggregate.
      const AggregateId& aggregate_id = aggregate_and_delta.first;

      const DemandAndFlowCount& demand_and_flow_count =
          nc::FindOrDieNoPrint(routing.demands(), aggregate_id);
      size_t flow_count = demand_and_flow_count.second;

      for (const FlowPathChange& change : aggregate_and_delta.second.changes) {
        uint64_t from_micros =
            duration_cast<microseconds>(change.from->delay()).count();
        uint64_t to_micros =
            duration_cast<microseconds>(change.to->delay()).count();

        int64_t delta = to_micros - from_micros;
        dist.Add(delta, flow_count);
      }
    }
  }

  absolute_path_stretch_micros->GetHandle(input.topology_file, input.tm_file,
                                          opt, new_link)
      ->AddValue(dist);
}

static void Record(const OptEvalInput& input,
                   const RoutingConfiguration& routing,
                   const RoutingConfiguration& new_routing,
                   const std::string& opt, bool new_link) {
  RoutingConfigurationDelta delta = routing.GetDifference(new_routing);
  double increasing_delay_count = 0;
  double decreasing_delay_count = 0;
  for (const auto& aggregate_and_delta : delta.aggregates) {
    const AggregateDelta& aggregate_delta = aggregate_and_delta.second;

    double on_longer_path = aggregate_delta.FractionOnLongerPath();
    double on_shorter_path = aggregate_delta.FractionDelta() -
                             aggregate_delta.FractionOnLongerPath();
    if (on_longer_path > 0.001) {
      ++increasing_delay_count;
    }

    if (on_shorter_path > 0.001) {
      ++decreasing_delay_count;
    }
  }

  double increasing_fraction = increasing_delay_count / delta.aggregates.size();
  double decreasing_fraction = decreasing_delay_count / delta.aggregates.size();
  increasing_delay_aggregate->GetHandle(input.topology_file, input.tm_file, opt,
                                        new_link)
      ->AddValue(increasing_fraction);
  decreasing_delay_aggreggate->GetHandle(input.topology_file, input.tm_file,
                                         opt, new_link)
      ->AddValue(decreasing_fraction);

  size_t overloaded_before = routing.OverloadedAggregates();
  size_t overloaded_after = routing.OverloadedAggregates();
  double net_overload = overloaded_after - overloaded_before;
  double overload_fraction = net_overload / delta.aggregates.size();
  increasing_overload_delta->GetHandle(input.topology_file, input.tm_file, opt,
                                       new_link)
      ->AddValue(overload_fraction);

  nc::net::Delay total_delay_before = routing.TotalPerFlowDelay();
  nc::net::Delay total_delay_after = new_routing.TotalPerFlowDelay();
  double change = (total_delay_after.count() - total_delay_before.count()) /
                  static_cast<double>(total_delay_before.count());
  total_delay_delta->GetHandle(input.topology_file, input.tm_file, opt,
                               new_link)
      ->AddValue(change);
}

static const nc::geo::CityData* GetCity(const std::string& name,
                                        nc::geo::Localizer* localizer) {
  nc::geo::FindCityRequest request;
  request.ascii_name = name;

  const nc::geo::CityData* city_data = localizer->FindCityOrNull(request);
  CHECK(city_data != nullptr);
  return city_data;
}

static void ProcessInput(const OptEvalInput& input, const std::string& opt,
                         nc::geo::Localizer* localizer) {
  const nc::lp::DemandMatrix& old_demand_matrix = *(input.demand_matrix);
  const nc::net::GraphStorage* old_graph = old_demand_matrix.graph();
  PathProvider old_path_provider(old_graph);
  auto old_routing = RunOptimizer(opt, old_demand_matrix, &old_path_provider);

  nc::net::GraphNodeSet all_nodes = old_graph->AllNodes();
  for (nc::net::GraphNodeIndex src : all_nodes) {
    for (nc::net::GraphNodeIndex dst : all_nodes) {
      if (src <= dst) {
        continue;
      }

      const std::string& src_id = old_graph->GetNode(src)->id();
      const std::string& dst_id = old_graph->GetNode(dst)->id();

      const nc::geo::CityData* src_city = GetCity(src_id, localizer);
      const nc::geo::CityData* dst_city = GetCity(dst_id, localizer);
      double distance_km = src_city->DistanceKm(*dst_city);
      uint32_t delay_micros = (distance_km / FLAGS_light_speed) * 1000.0;
      delay_micros = std::max(delay_micros, static_cast<uint32_t>(1));
      bool new_link = !old_graph->HasLink(src_id, dst_id);

      nc::net::Bandwidth bw =
          nc::net::Bandwidth::FromMBitsPerSecond(FLAGS_link_speed_Mbps);

      nc::net::GraphBuilder new_topology = old_graph->ToBuilder();
      new_topology.AddLink({src_id, dst_id, bw, microseconds(delay_micros)});
      new_topology.AddLink({dst_id, src_id, bw, microseconds(delay_micros)});
      new_topology.RemoveMultipleLinks();
      nc::net::GraphStorage new_graph(new_topology);
      PathProvider new_path_provider(&new_graph);
      auto new_routing =
          RunOptimizer(opt, nc::lp::DemandMatrix(old_demand_matrix, &new_graph),
                       &new_path_provider);
      Record(input, *old_routing, *new_routing, opt, new_link);
      RecordStretch(input, *old_routing, {new_routing.get()}, opt, new_link);
    }
  }
}

}  // namespace ctr

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::metrics::InitMetrics();
  nc::geo::Localizer localizer(FLAGS_cities_file);

  // Will switch off timestamps.
  auto timestamp_provider =
      ::nc::make_unique<nc::metrics::NullTimestampProvider>();
  nc::metrics::DefaultMetricManager()->set_timestamp_provider(
      std::move(timestamp_provider));

  std::vector<std::unique_ptr<nc::net::GraphStorage>> graphs;
  std::vector<ctr::OptEvalInput> to_process;
  std::tie(graphs, to_process) = ctr::GetOptEvalInputs();

  std::vector<std::string> optimizers = nc::Split(FLAGS_optimizers, ",");
  for (const auto& opt : optimizers) {
    for (const auto& input : to_process) {
      ctr::ProcessInput(input, opt, &localizer);
    }
  }
}
