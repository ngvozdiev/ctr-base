#ifndef CTR_OPT_COMPARE_H
#define CTR_OPT_COMPARE_H

#include "../common.h"
#include "opt.h"

namespace ctr {

// Performance of a single optimizer on a single TM.
class PerOptimizerSummary {
 public:
  PerOptimizerSummary(const TrafficMatrix& input,
                      const RoutingConfiguration& output,
                      uint64_t processing_time_ms,
                      const nc::net::GraphStorage* graph)
      : processing_time_ms(processing_time_ms) {
    // How far away each flow is from the shortest path. In seconds.
    std::vector<double> path_stretches;
    std::vector<double> path_stretches_rel;

    // For each link the total flow over the link and the link's capacity.
    std::vector<std::pair<double, double>> link_load;

    // Number of paths per aggregate.
    std::vector<double> num_paths;

    // A map from a link to the total load over the link.
    std::map<nc::net::GraphLinkIndex, uint64_t> link_to_total_load;

    // Unmet demand.
    std::vector<double> unmet_demand;

    for (const auto& aggregate_and_aggregate_output : output.routes()) {
      const AggregateId& aggregate_id = aggregate_and_aggregate_output.first;
      const std::vector<RouteAndFraction>& routes =
          aggregate_and_aggregate_output.second;
      const DemandAndFlowCount& demand_and_flow_count =
          nc::FindOrDie(input.demands(), aggregate_id);

      size_t total_num_flows = demand_and_flow_count.second;
      double total_aggregate_demand = demand_and_flow_count.first.bps();
      double required_per_flow_bps = total_aggregate_demand / total_num_flows;

      std::unique_ptr<nc::net::Walk> sp = nc::net::ShortestPathWithConstraints(
          aggregate_id.src(), aggregate_id.dst(), *graph, {});
      CHECK(sp);

      std::chrono::microseconds shortest_path_delay = sp->delay();
      double sp_delay_sec =
          std::chrono::duration<double>(shortest_path_delay).count();

      // Will limit the delay at 1ms.
      sp_delay_sec = std::max(sp_delay_sec, 0.001);

      size_t path_count = 0;
      for (const auto& route : routes) {
        const nc::net::Walk* path = route.first;
        double fraction = route.second;

        const fubar::PathOutput& path_output = tag_and_path.second;
        std::chrono::microseconds path_delay = path->delay();
        CHECK(path_delay >= shortest_path_delay);

        if (path_output.bps_limit() > 0) {
          ++path_count;
        }

        double path_delay_sec =
            std::chrono::duration<double>(path_delay).count();
        path_delay_sec = std::max(path_delay_sec, 0.001);
        double delay_sec = path_delay_sec - sp_delay_sec;
        double unmet = std::max(0.0, required_per_flow_bps - per_flow_bps);
        double delta_rel = path_delay_sec / sp_delay_sec;

        for (size_t i = 0; i < path_output.fraction() * total_num_flows; ++i) {
          path_stretches.emplace_back(delay_sec);
          path_stretches_rel.emplace_back(delta_rel);
          unmet_demand.emplace_back(unmet);
        }

        for (ncode::net::GraphLinkIndex link :
             path_output.path()->link_sequence().links()) {
          link_to_total_load[link] +=
              total_aggregate_demand * path_output.fraction();
        }
      }

      num_paths.emplace_back(path_count);
    }

    ncode::net::GraphLinkSet links_seen;
    for (const auto& link_and_total_load : link_to_total_load) {
      ncode::net::GraphLinkIndex link_index = link_and_total_load.first;
      links_seen.Insert(link_index);
      const ncode::net::GraphLink* link =
          path_cache->graph_storage()->GetLink(link_index);

      double total_load = link_and_total_load.second;
      link_load.emplace_back(total_load, link->bandwidth().bps());
    }

    for (ncode::net::GraphLinkIndex link_index :
         path_cache->graph_storage()->AllLinks()) {
      if (!links_seen.Contains(link_index)) {
        const ncode::net::GraphLink* link =
            path_cache->graph_storage()->GetLink(link_index);
        link_load.emplace_back(0, link->bandwidth().bps());
      }
    }
  }

  void SaveStats(const std::string& file) const {
    std::string prefix = nc::StrCat(aggregate_count, ",");

    std::string all_stretches = nc::Join(path_stretches, ",");
    nc::File::WriteStringToFileOrDie(nc::StrCat(prefix, all_stretches),
                                     nc::StrCat(file, "_stretches"));

    std::string all_stretches_rel = nc::Join(path_stretches_rel, ",");
    nc::File::WriteStringToFileOrDie(nc::StrCat(prefix, all_stretches_rel),
                                     nc::StrCat(file, "_stretches_rel"));

    std::function<std::string(
        const std::pair<double, double>& link_flow_and_capacity)> f =
        [](const std::pair<double, double>& link_flow_and_capacity) {
          return nc::StrCat(link_flow_and_capacity.first, ",",
                            link_flow_and_capacity.second);
        };
    std::string all_link_utilization = nc::Join(link_loads, ",", f);
    nc::File::WriteStringToFileOrDie(nc::StrCat(prefix, all_link_utilization),
                                     nc::StrCat(file, "_link_utilization"));

    std::string all_unmet_demand = nc::Join(unmet_demand, ",");
    nc::File::WriteStringToFileOrDie(nc::StrCat(prefix, all_unmet_demand),
                                     nc::StrCat(file, "_unmet_demand"));

    std::string all_num_paths = nc::Join(num_paths, ",");
    nc::File::WriteStringToFileOrDie(nc::StrCat(prefix, all_num_paths),
                                     nc::StrCat(file, "_num_paths"));
  }

  void SaveProcessingTime(const std::string& file) const {
    std::string prefix = nc::StrCat(aggregate_count, ",");

    nc::File::WriteStringToFileOrDie(
        nc::StrCat(prefix, std::to_string(processing_time_ms)),
        nc::StrCat(file, "_processing_time"));

    if (network_state_path_count != std::numeric_limits<uint64_t>::max()) {
      nc::File::WriteStringToFileOrDie(
          std::to_string(network_state_path_count),
          nc::StrCat(file, "_network_state_path_count"));
    }
  }

  std::vector<double> LinkUtilization() {
    std::vector<double> out;
    for (const std::pair<double, double>& link_flow_and_capacity : link_loads) {
      out.emplace_back(link_flow_and_capacity.first /
                       link_flow_and_capacity.second);
    }

    return out;
  }

 private:
  // Per-flow path stretch.
  std::vector<double> path_stretches;
  std::vector<double> path_stretches_rel;

  // Link utilization, first is flow over link, second is link capacity.
  std::vector<std::pair<double, double>> link_loads;

  // Number of paths per aggregate.
  std::vector<double> num_paths;

  // Unmet demand
  std::vector<double> unmet_demand;

  // Time to run the optimizer.
  uint64_t processing_time_ms;

  // The number of links that are oversubscribed.
  uint64_t oversubscribed_links;

  // Total number of paths that need to be installed in the network, this
  // includes alternative paths to protect against failures, as well as primary
  // paths. Paths that are the same as the shortest path are not included.
  size_t network_state_path_count;

  // Number of aggregates.
  uint64_t aggregate_count;
};

}  // namespace ctr

#endif
