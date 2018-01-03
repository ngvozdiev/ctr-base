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

#include "ncode/common.h"
#include "ncode/file.h"
#include "ncode/logging.h"
#include "ncode/lp/demand_matrix.h"
#include "ncode/net/net_common.h"
#include "ncode/strutil.h"
#include "ncode/viz/web_page.h"
#include "ncode/thread_runner.h"
#include "common.h"
#include "geo/geo.h"
#include "metrics/metrics.h"
#include "opt/ctr.h"
#include "opt/opt.h"
#include "opt/path_provider.h"
#include "demand_matrix_input.h"

DEFINE_string(cities_file, "cities5000.txt", "Location of the cities file");
DEFINE_double(
    light_speed, 200.0,
    "Speed of light (in km per millisecond), default is speed in fiber");
DEFINE_double(link_speed_Mbps, 1000, "All new links will be this fast");
DEFINE_string(optimizers, "SP,CTR,CTRNFC,MinMax,MinMaxLD,MinMaxK10,B4",
              "The optimizers to use, comma-separated.");
DEFINE_bool(fixed_flow_count, false,
            "If true all aggregates will have the same flow count, if false "
            "the flow count will depend on the size of the aggregate");
DEFINE_uint64(threads, 4, "Number of parallel threads to run");

using namespace std::chrono;

static auto* decreasing_delay_aggreggate =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string,
                               bool, std::string>(
            "decreasing_delay_aggregate_fraction",
            "Fraction of aggregates whose delay decreases", "Topology",
            "Traffic matrix", "Optimizer", "New link",
            "Relative stretch threshold");

static auto* increasing_delay_aggregate =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string,
                               bool, std::string>(
            "increasing_delay_aggregate_fraction",
            "Fraction of aggregates whose delay increases", "Topology",
            "Traffic matrix", "Optimizer", "New link",
            "Relative stretch threshold");

static auto* decreasing_delay_flows =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string,
                               bool, std::string>(
            "decreasing_delay_flows_fraction",
            "Fraction of flows whose delay decreases", "Topology",
            "Traffic matrix", "Optimizer", "New link",
            "Relative stretch threshold");

static auto* increasing_delay_flows =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string,
                               bool, std::string>(
            "increasing_delay_flows_fraction",
            "Fraction of flows whose delay increases", "Topology",
            "Traffic matrix", "Optimizer", "New link",
            "Relative stretch threshold");

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

static auto* total_delay =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<double, std::string, std::string, std::string,
                               bool>("total_absolute_delay", "Total delay",
                                     "Topology", "Traffic matrix", "Optimizer",
                                     "New link");

static auto* absolute_path_stretch_ms =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<nc::DiscreteDistribution<int64_t>, std::string,
                               std::string, std::string, bool>(
            "absolute_path_stretch_ms",
            "Distribution of absolute per-flow deltas in milliseconds",
            "Topology", "Traffic matrix", "Optimizer", "New link");

static auto* relative_path_stretch =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<nc::DiscreteDistribution<int64_t>, std::string,
                               std::string, std::string, bool>(
            "relative_path_stretch",
            "Distribution of relative path stretch (quantized x10000)",
            "Topology", "Traffic matrix", "Optimizer", "New link");

static auto* link_delay_micros =
    nc::metrics::DefaultMetricManager()
        -> GetThreadSafeMetric<uint64_t, std::string, std::string, std::string,
                               bool>("link_delay_micros",
                                     "Link delay in microseconds", "Topology",
                                     "Traffic matrix", "Optimizer", "New link");

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

    if (FLAGS_fixed_flow_count) {
      flow_count = kTopAggregateFlowCount;
    }

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
    const DemandMatrixAndFilename& input, const RoutingConfiguration& routing,
    const std::vector<const RoutingConfiguration*> new_routings,
    const std::string& opt, bool new_link) {
  nc::DiscreteDistribution<int64_t> dist;
  nc::DiscreteDistribution<int64_t> dist_relative;

  for (const RoutingConfiguration* new_routing : new_routings) {
    RoutingConfigurationDelta delta = routing.GetDifference(*new_routing);
    for (const auto& aggregate_and_delta : delta.aggregates) {
      // Need to figure out how many flows there are in an aggregate.
      const AggregateId& aggregate_id = aggregate_and_delta.first;

      const DemandAndFlowCount& demand_and_flow_count =
          nc::FindOrDieNoPrint(routing.demands(), aggregate_id);
      size_t flow_count = demand_and_flow_count.second;

      for (const FlowPathChange& change : aggregate_and_delta.second.changes) {
        int64_t from_ms =
            duration_cast<milliseconds>(change.from->delay()).count();
        int64_t to_ms = duration_cast<milliseconds>(change.to->delay()).count();
        size_t flows_changed = change.fraction * flow_count;
        if (flows_changed == 0) {
          continue;
        }

        from_ms = std::max(from_ms, static_cast<int64_t>(1));
        to_ms = std::max(to_ms, static_cast<int64_t>(1));

        int64_t delta = to_ms - from_ms;
        dist.Add(delta, flows_changed);

        double delta_rel = delta / static_cast<double>(from_ms);
        dist_relative.Add(delta_rel * 10000, flows_changed);
      }
    }
  }

  absolute_path_stretch_ms->GetHandle(input.topology_file, input.file, opt,
                                      new_link)
      ->AddValue(dist);
  relative_path_stretch->GetHandle(input.topology_file, input.file, opt,
                                   new_link)
      ->AddValue(dist_relative);
}

static void Record(const DemandMatrixAndFilename& input,
                   const RoutingConfiguration& routing,
                   const RoutingConfiguration& new_routing,
                   const std::string& opt, bool new_link, double threshold) {
  RoutingConfigurationDelta delta = routing.GetDifference(new_routing);
  double increasing_delay_count = 0;
  double decreasing_delay_count = 0;
  double increasing_flow_count = 0;
  double decreasing_flow_count = 0;
  double total_flow_count = 0;

  for (const auto& aggregate_and_delta : delta.aggregates) {
    const AggregateDelta& aggregate_delta = aggregate_and_delta.second;

    double on_longer_path = aggregate_delta.FractionOnLongerPath(threshold);
    double on_shorter_path =
        aggregate_delta.FractionDelta(threshold) - on_longer_path;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(routing.demands(), aggregate_and_delta.first);
    double flow_count = demand_and_flow_count.second;
    total_flow_count += flow_count;
    increasing_flow_count += on_longer_path * flow_count;
    decreasing_flow_count += on_shorter_path * flow_count;

    if (on_longer_path > 0.001) {
      ++increasing_delay_count;
    }

    if (on_shorter_path > 0.001) {
      ++decreasing_delay_count;
    }
  }

  std::string threshold_str = nc::ToStringMaxDecimals(threshold, 1);
  double increasing_flow_fraction = increasing_flow_count / total_flow_count;
  double decreasing_flow_fraction = decreasing_flow_count / total_flow_count;
  increasing_delay_flows->GetHandle(input.topology_file, input.file, opt,
                                    new_link, threshold_str)
      ->AddValue(increasing_flow_fraction);
  decreasing_delay_flows->GetHandle(input.topology_file, input.file, opt,
                                    new_link, threshold_str)
      ->AddValue(decreasing_flow_fraction);

  double increasing_fraction = increasing_delay_count / delta.aggregates.size();
  double decreasing_fraction = decreasing_delay_count / delta.aggregates.size();
  increasing_delay_aggregate->GetHandle(input.topology_file, input.file, opt,
                                        new_link, threshold_str)
      ->AddValue(increasing_fraction);
  decreasing_delay_aggreggate->GetHandle(input.topology_file, input.file, opt,
                                         new_link, threshold_str)
      ->AddValue(decreasing_fraction);
}

static void RecordOverloadAndTotal(const DemandMatrixAndFilename& input,
                                   const RoutingConfiguration& routing,
                                   const RoutingConfiguration& new_routing,
                                   const std::string& opt, bool new_link) {
  size_t overloaded_before = routing.OverloadedAggregates();
  size_t overloaded_after = routing.OverloadedAggregates();
  double net_overload = overloaded_after - overloaded_before;
  double overload_fraction = net_overload / routing.demands().size();
  increasing_overload_delta->GetHandle(input.topology_file, input.file, opt,
                                       new_link)
      ->AddValue(overload_fraction);

  nc::net::Delay total_delay_before = routing.TotalPerFlowDelay();
  nc::net::Delay total_delay_after = new_routing.TotalPerFlowDelay();
  double change = (total_delay_after.count() - total_delay_before.count()) /
                  static_cast<double>(total_delay_before.count());
  total_delay_delta->GetHandle(input.topology_file, input.file, opt, new_link)
      ->AddValue(change);
  total_delay->GetHandle(input.topology_file, input.file, opt, new_link)
      ->AddValue(total_delay_after.count());
}

static void RecordLinkDelay(const DemandMatrixAndFilename& input,
                            std::chrono::microseconds delay,
                            const std::string& opt, bool new_link) {
  link_delay_micros->GetHandle(input.topology_file, input.file, opt, new_link)
      ->AddValue(delay.count());
}

static const nc::geo::CityData* GetCity(const std::string& name,
                                        nc::geo::Localizer* localizer) {
  nc::geo::FindCityRequest request;
  request.ascii_name = name;

  const nc::geo::CityData* city_data = localizer->FindCityOrNull(request);
  CHECK(city_data != nullptr);
  return city_data;
}

static void ProcessInput(const DemandMatrixAndFilename& input,
                         const std::string& opt,
                         nc::geo::Localizer* localizer) {
  LOG(INFO) << "Processing " << input.topology_file << " tm " << input.file
            << " opt " << opt;
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
      for (double threshold = 0.0; threshold <= 1.0; threshold += 0.1) {
        Record(input, *old_routing, *new_routing, opt, new_link, threshold);
      }
      RecordOverloadAndTotal(input, *old_routing, *new_routing, opt, new_link);
      RecordStretch(input, *old_routing, {new_routing.get()}, opt, new_link);
      RecordLinkDelay(input, microseconds(delay_micros), opt, new_link);
    }
  }
  LOG(INFO) << "Done";
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

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> to_process;
  std::tie(topologies, to_process) = ctr::GetDemandMatrixInputs();

  std::vector<std::string> optimizers = nc::Split(FLAGS_optimizers, ",");
  for (const auto& opt : optimizers) {
    nc::RunInParallel<ctr::DemandMatrixAndFilename>(
        to_process,
        [&opt, &localizer](const ctr::DemandMatrixAndFilename& input) {
          ctr::ProcessInput(input, opt, &localizer);
        },
        FLAGS_threads);

    //    for (const auto& input : to_process) {
    //      ctr::ProcessInput(input, opt, &localizer);
    //    }
  }
}
