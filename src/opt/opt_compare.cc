#include "opt_compare.h"

namespace ctr {

std::string DumpRoutingInfo(const RoutingConfiguration& routing) {
  const nc::net::GraphStorage* graph = routing.graph();

  // A map from a link to the total load over the link.
  std::map<nc::net::GraphLinkIndex, nc::net::Bandwidth> link_to_total_load;

  OverSubModel model(routing);
  const std::map<const nc::net::Walk*, nc::net::Bandwidth>& per_flow_rates =
      model.per_flow_bandwidth_map();

  std::vector<double> path_stretches;
  std::vector<double> path_stretches_rel;
  std::vector<double> unmet_demand;
  std::vector<uint32_t> num_paths;

  for (const auto& aggregate_and_aggregate_output : routing.routes()) {
    const AggregateId& aggregate_id = aggregate_and_aggregate_output.first;
    const std::vector<RouteAndFraction>& routes =
        aggregate_and_aggregate_output.second;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDie(routing.demands(), aggregate_id);

    size_t total_num_flows = demand_and_flow_count.second;
    nc::net::Bandwidth total_aggregate_demand = demand_and_flow_count.first;
    nc::net::Bandwidth required_per_flow =
        total_aggregate_demand / total_num_flows;

    std::chrono::microseconds shortest_path_delay =
        aggregate_id.GetSPDelay(*graph);
    double sp_delay_sec =
        std::chrono::duration<double>(shortest_path_delay).count();

    // Will limit the delay at 1ms.
    sp_delay_sec = std::max(sp_delay_sec, 0.001);

    size_t path_count = 0;
    for (const auto& route : routes) {
      const nc::net::Walk* path = route.first;
      double fraction = route.second;
      CHECK(fraction > 0);

      std::chrono::microseconds path_delay = path->delay();
      CHECK(path_delay >= shortest_path_delay);

      ++path_count;

      double path_delay_sec = std::chrono::duration<double>(path_delay).count();
      path_delay_sec = std::max(path_delay_sec, 0.001);
      double delay_sec = path_delay_sec - sp_delay_sec;

      nc::net::Bandwidth per_flow_rate = per_flow_rates[path];
      double unmet =
          std::max(0.0, required_per_flow.Mbps() - per_flow_rate.Mbps());
      double delta_rel = path_delay_sec / sp_delay_sec;

      for (size_t i = 0; i < fraction * total_num_flows; ++i) {
        path_stretches.emplace_back(delay_sec);
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

    double total_load = link_and_total_load.second;
    link_loads.emplace_back(total_load / link->bandwidth());
  }

  for (nc::net::GraphLinkIndex link_index : graph->AllLinks()) {
    if (!links_seen.Contains(link_index)) {
      link_loads.emplace_back(0);
    }
  }

  std::string out;
  nc::StrAppend(&out, std::to_string(routing.routes().size()), " ");
  nc::StrAppend(&out, nc::Join(path_stretches, " "), " ");
  nc::StrAppend(&out, nc::Join(path_stretches_rel, " "), " ");
  nc::StrAppend(&out, nc::Join(num_paths, " "), " ");
}

}  // namespace ctr
