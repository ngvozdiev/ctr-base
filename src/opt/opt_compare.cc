#include "opt_compare.h"

#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/perfect_hash.h"
#include "ncode_common/src/stats.h"
#include "ncode_common/src/strutil.h"
#include "../common.h"
#include "oversubscription_model.h"

namespace ctr {

static nc::viz::NpyArray::Types GetRCTypes() {
  using namespace nc::viz;
  NpyArray::Types types = {{"num_aggregates", NpyArray::UINT32}};
  for (uint32_t i = 0; i < 101; ++i) {
    types.emplace_back(nc::StrCat("path_stretch_ms_p", i), NpyArray::UINT32);
  }
  for (uint32_t i = 0; i < 101; ++i) {
    types.emplace_back(nc::StrCat("path_stretch_rel_p", i), NpyArray::UINT32);
  }
  for (uint32_t i = 0; i < 101; ++i) {
    types.emplace_back(nc::StrCat("paths_per_aggregate_p", i),
                       NpyArray::UINT32);
  }
  for (uint32_t i = 0; i < 101; ++i) {
    types.emplace_back(nc::StrCat("link_utilization_p", i), NpyArray::DOUBLE);
  }
  for (uint32_t i = 0; i < 101; ++i) {
    types.emplace_back(nc::StrCat("unmet_demand_mbps_p", i), NpyArray::DOUBLE);
  }

  return types;
}

RoutingConfigInfo::RoutingConfigInfo() : nc::viz::NpyArray(GetRCTypes()) {}

void RoutingConfigInfo::Add(const RoutingConfiguration& routing) {
  using namespace std::chrono;
  const nc::net::GraphStorage* graph = routing.graph();

  // A map from a link to the total load over the link.
  std::map<nc::net::GraphLinkIndex, nc::net::Bandwidth> link_to_total_load;

  OverSubModel model(routing);
  const std::map<const nc::net::Walk*, nc::net::Bandwidth>& per_flow_rates =
      model.per_flow_bandwidth_map();

  std::vector<milliseconds> path_stretches;
  std::vector<double> path_stretches_rel;
  std::vector<nc::net::Bandwidth> unmet_demand;
  std::vector<size_t> num_paths;

  for (const auto& aggregate_and_aggregate_output : routing.routes()) {
    const AggregateId& aggregate_id = aggregate_and_aggregate_output.first;
    const std::vector<RouteAndFraction>& routes =
        aggregate_and_aggregate_output.second;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(routing.demands(), aggregate_id);

    size_t total_num_flows = demand_and_flow_count.second;
    nc::net::Bandwidth total_aggregate_demand = demand_and_flow_count.first;
    nc::net::Bandwidth required_per_flow =
        total_aggregate_demand / total_num_flows;

    microseconds shortest_path_delay = aggregate_id.GetSPDelay(*graph);
    milliseconds sp_delay_ms = duration_cast<milliseconds>(shortest_path_delay);

    // Will limit the delay at 1ms.
    sp_delay_ms = std::max(sp_delay_ms, milliseconds(1));

    size_t path_count = 0;
    for (const auto& route : routes) {
      const nc::net::Walk* path = route.first;
      double fraction = route.second;
      CHECK(fraction > 0);

      microseconds path_delay = path->delay();
      CHECK(path_delay >= shortest_path_delay);

      ++path_count;

      milliseconds path_delay_ms = duration_cast<milliseconds>(path_delay);
      path_delay_ms = std::max(path_delay_ms, milliseconds(1));
      milliseconds delay_delta = path_delay_ms - sp_delay_ms;

      nc::net::Bandwidth per_flow_rate =
          nc::FindOrDieNoPrint(per_flow_rates, path);
      nc::net::Bandwidth unmet = std::max(nc::net::Bandwidth::Zero(),
                                          required_per_flow - per_flow_rate);
      double delta_rel =
          static_cast<double>(path_delay_ms.count()) / sp_delay_ms.count();

      for (size_t i = 0; i < fraction * total_num_flows; ++i) {
        path_stretches.emplace_back(delay_delta);
        path_stretches_rel.emplace_back(delta_rel);
        unmet_demand.emplace_back(unmet);
      }

      for (nc::net::GraphLinkIndex link : path->links()) {
        link_to_total_load[link] += total_aggregate_demand * fraction;
      }
    }

    num_paths.emplace_back(path_count);
  }

  std::vector<double> link_loads;
  nc::net::GraphLinkSet links_seen;
  for (const auto& link_and_total_load : link_to_total_load) {
    nc::net::GraphLinkIndex link_index = link_and_total_load.first;
    links_seen.Insert(link_index);
    const nc::net::GraphLink* link = graph->GetLink(link_index);

    nc::net::Bandwidth total_load = link_and_total_load.second;
    link_loads.emplace_back(total_load / link->bandwidth());
  }

  for (nc::net::GraphLinkIndex link_index : graph->AllLinks()) {
    if (!links_seen.Contains(link_index)) {
      link_loads.emplace_back(0);
    }
  }

  std::vector<milliseconds> path_stretches_p = nc::Percentiles(&path_stretches);
  std::vector<double> path_stretches_rel_p =
      nc::Percentiles(&path_stretches_rel);
  std::vector<nc::net::Bandwidth> unmet_demand_p =
      nc::Percentiles(&unmet_demand);
  std::vector<size_t> num_paths_p = nc::Percentiles(&num_paths);
  std::vector<double> link_loads_p = nc::Percentiles(&link_loads);

  std::vector<nc::viz::NpyArray::StringOrNumeric> row = {
      routing.routes().size()};
  for (milliseconds v : path_stretches_p) {
    row.emplace_back(v.count());
  }
  for (double v : path_stretches_rel_p) {
    row.emplace_back(v);
  }
  for (size_t v : num_paths_p) {
    row.emplace_back(v);
  }
  for (double v : link_loads_p) {
    row.emplace_back(v);
  }
  for (nc::net::Bandwidth v : unmet_demand_p) {
    row.emplace_back(v.Mbps());
  }

  AddRow(row);
}

static nc::viz::NpyArray::Types GetRCDTypes() {
  using namespace nc::viz;
  NpyArray::Types types = {{"volume_delta", NpyArray::DOUBLE},
                           {"flow_count_delta", NpyArray::DOUBLE},
                           {"volume_delta_longer_path", NpyArray::DOUBLE},
                           {"flow_count_delta_longer_path", NpyArray::DOUBLE},
                           {"add_count", NpyArray::UINT32},
                           {"update_count", NpyArray::UINT32},
                           {"remove_count", NpyArray::UINT32}};
  for (uint32_t i = 0; i < 101; ++i) {
    types.emplace_back(nc::StrCat("per_aggregate_fraction_p", i),
                       NpyArray::DOUBLE);
  }

  return types;
}

RoutingConfigDeltaInfo::RoutingConfigDeltaInfo()
    : nc::viz::NpyArray(GetRCDTypes()) {}

void RoutingConfigDeltaInfo::Add(const RoutingConfigurationDelta& delta) {
  using namespace std::chrono;

  std::vector<double> fraction_deltas;
  for (const auto& aggregate_id_and_delta : delta.aggregates) {
    const ctr::AggregateDelta& aggregate_delta = aggregate_id_and_delta.second;
    fraction_deltas.emplace_back(aggregate_delta.fraction_delta);
  }
  std::vector<double> fraction_deltas_p = nc::Percentiles(&fraction_deltas);

  size_t route_adds;
  size_t route_removals;
  size_t route_updates;
  std::tie(route_adds, route_removals, route_updates) = delta.TotalRoutes();

  std::vector<nc::viz::NpyArray::StringOrNumeric> row = {
      delta.total_volume_fraction_delta,
      delta.total_flow_fraction_delta,
      delta.total_volume_fraction_on_longer_path,
      delta.total_flow_fraction_on_longer_path,
      route_adds,
      route_updates,
      route_removals};
  for (double v : fraction_deltas_p) {
    row.emplace_back(v);
  }

  AddRow(row);
}

}  // namespace ctr
