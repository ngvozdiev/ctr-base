#include "opt_compare.h"

namespace ctr {

PerOptimizerSummary::PerOptimizerSummary(const TrafficMatrix& input,
                                         const RoutingConfiguration& output,
                                         double link_scale_factor,
                                         uint64_t processing_time_ms)
    : processing_time_ms_(processing_time_ms),
      aggregate_count_(input.demands().size()) {
  const nc::net::GraphStorage* graph = input.graph();

  // A map from a link to the total load over the link.
  std::map<nc::net::GraphLinkIndex, uint64_t> link_to_total_load;

  OverSubModel model(input, output, link_scale_factor);
  const std::map<const nc::net::Walk*, nc::net::Bandwidth>& per_flow_rates =
      model.per_flow_bandwidth_map();

  for (const auto& aggregate_and_aggregate_output : output.routes()) {
    const AggregateId& aggregate_id = aggregate_and_aggregate_output.first;
    const std::vector<RouteAndFraction>& routes =
        aggregate_and_aggregate_output.second;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDie(input.demands(), aggregate_id);

    size_t total_num_flows = demand_and_flow_count.second;
    double total_aggregate_demand = demand_and_flow_count.first.bps();
    double required_per_flow_bps = total_aggregate_demand / total_num_flows;

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
      double unmet = std::max(0.0, required_per_flow_bps - per_flow_rate.bps());
      double delta_rel = path_delay_sec / sp_delay_sec;

      for (size_t i = 0; i < fraction * total_num_flows; ++i) {
        path_stretches_.emplace_back(delay_sec);
        path_stretches_rel_.emplace_back(delta_rel);
        unmet_demand_.emplace_back(unmet);
      }

      for (nc::net::GraphLinkIndex link : path->links()) {
        link_to_total_load[link] += total_aggregate_demand * fraction;
      }
    }

    num_paths_.emplace_back(path_count);
  }

  nc::net::GraphLinkSet links_seen;
  for (const auto& link_and_total_load : link_to_total_load) {
    nc::net::GraphLinkIndex link_index = link_and_total_load.first;
    links_seen.Insert(link_index);
    const nc::net::GraphLink* link = graph->GetLink(link_index);

    double total_load = link_and_total_load.second;
    link_loads_.emplace_back(total_load, link->bandwidth().bps());
  }

  for (nc::net::GraphLinkIndex link_index : graph->AllLinks()) {
    if (!links_seen.Contains(link_index)) {
      const nc::net::GraphLink* link = graph->GetLink(link_index);
      link_loads_.emplace_back(0, link->bandwidth().bps());
    }
  }
}

void PerOptimizerSummary::SaveStats(const std::string& file) const {
  std::string prefix = nc::StrCat(aggregate_count_, ",");

  std::string all_stretches = nc::Join(path_stretches_, ",");
  nc::File::WriteStringToFileOrDie(nc::StrCat(prefix, all_stretches),
                                   nc::StrCat(file, "_stretches"));

  std::string all_stretches_rel = nc::Join(path_stretches_rel_, ",");
  nc::File::WriteStringToFileOrDie(nc::StrCat(prefix, all_stretches_rel),
                                   nc::StrCat(file, "_stretches_rel"));

  std::function<std::string(
      const std::pair<double, double>& link_flow_and_capacity)> f =
      [](const std::pair<double, double>& link_flow_and_capacity) {
        return nc::StrCat(link_flow_and_capacity.first, ",",
                          link_flow_and_capacity.second);
      };
  std::string all_link_utilization = nc::Join(link_loads_, ",", f);
  nc::File::WriteStringToFileOrDie(nc::StrCat(prefix, all_link_utilization),
                                   nc::StrCat(file, "_link_utilization"));

  std::string all_unmet_demand = nc::Join(unmet_demand_, ",");
  nc::File::WriteStringToFileOrDie(nc::StrCat(prefix, all_unmet_demand),
                                   nc::StrCat(file, "_unmet_demand"));

  std::string all_num_paths = nc::Join(num_paths_, ",");
  nc::File::WriteStringToFileOrDie(nc::StrCat(prefix, all_num_paths),
                                   nc::StrCat(file, "_num_paths"));

  nc::File::WriteStringToFileOrDie(
      nc::StrCat(prefix, std::to_string(processing_time_ms_)),
      nc::StrCat(file, "_processing_time"));
}

std::vector<double> PerOptimizerSummary::LinkUtilization() const {
  std::vector<double> out;
  for (const std::pair<double, double>& link_flow_and_capacity : link_loads_) {
    out.emplace_back(link_flow_and_capacity.first /
                     link_flow_and_capacity.second);
  }

  return out;
}

}  // namespace ctr
