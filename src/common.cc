#include "ncode_common/src/common.h"
#include "common.h"

#include <tuple>

#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/net/algorithm.h"
#include "ncode_common/src/viz/graph.h"
#include "opt/oversubscription_model.h"

namespace ctr {

bool operator<(const AggregateId& a, const AggregateId& b) {
  return std::tie(a.src_, a.dst_) < std::tie(b.src_, b.dst_);
}

bool operator==(const AggregateId& a, const AggregateId& b) {
  return std::tie(a.src_, a.dst_) == std::tie(b.src_, b.dst_);
}

bool operator!=(const AggregateId& a, const AggregateId& b) {
  return std::tie(a.src_, a.dst_) != std::tie(b.src_, b.dst_);
}

std::string AggregateId::ToString(const nc::net::GraphStorage& graph) const {
  std::string src_str = graph.GetNode(src_)->id();
  std::string dst_str = graph.GetNode(dst_)->id();
  return nc::StrCat("<", src_str, ",", dst_str, ">");
}

nc::net::Delay AggregateId::GetSPDelay(
    const nc::net::GraphStorage& graph) const {
  std::unique_ptr<nc::net::Walk> sp =
      nc::net::ShortestPathWithConstraints(src_, dst_, graph, {});
  CHECK(sp);
  return sp->delay();
}

void TrafficMatrix::AddDemand(const AggregateId& aggregate_id,
                              const DemandAndFlowCount& demand_and_flow_count) {
  CHECK(!nc::ContainsKey(demands_, aggregate_id));
  demands_[aggregate_id] = demand_and_flow_count;
}

TrafficMatrix::TrafficMatrix(
    const nc::lp::DemandMatrix& demand_matrix,
    const std::map<nc::lp::SrcAndDst, size_t>& flow_counts)
    : graph_(demand_matrix.graph()) {
  for (const auto& element : demand_matrix.elements()) {
    size_t flow_count =
        nc::FindWithDefault(flow_counts, {element.src, element.dst}, 1ul);
    demands_[{element.src, element.dst}] = {element.demand, flow_count};
  }
}

std::unique_ptr<nc::lp::DemandMatrix> TrafficMatrix::ToDemandMatrix() const {
  std::vector<nc::lp::DemandMatrixElement> elements;
  for (const auto& aggregate_and_demand : demands_) {
    const AggregateId& aggregate = aggregate_and_demand.first;
    nc::net::Bandwidth demand = aggregate_and_demand.second.first;

    elements.emplace_back(aggregate.src(), aggregate.dst(), demand);
  }

  return nc::make_unique<nc::lp::DemandMatrix>(std::move(elements), graph_);
}

std::unique_ptr<TrafficMatrix> TrafficMatrix::ScaleDemands(
    double factor, const std::set<AggregateId>& to_scale) const {
  std::map<AggregateId, DemandAndFlowCount> new_demands;
  for (const auto& aggregate_and_demand : demands_) {
    const AggregateId& aggregate = aggregate_and_demand.first;
    const DemandAndFlowCount& demand_and_flow_count =
        aggregate_and_demand.second;
    if (!to_scale.empty() && !nc::ContainsKey(to_scale, aggregate)) {
      new_demands[aggregate] = demand_and_flow_count;
      continue;
    }

    new_demands[aggregate] = {demand_and_flow_count.first * factor,
                              demand_and_flow_count.second};
  }

  return nc::make_unique<TrafficMatrix>(graph_, new_demands);
}

static double Pick(double current, double fraction, std::mt19937* rnd) {
  double delta = current * fraction;
  double low = current - delta;
  double high = current + delta;
  std::uniform_real_distribution<double> dist(low, high);
  double new_demand = dist(*rnd);
  return std::max(new_demand, 1.0);
}

std::unique_ptr<TrafficMatrix> TrafficMatrix::Randomize(
    double demand_fraction, double flow_count_fraction, size_t aggregate_count,
    std::mt19937* rnd) const {
  std::vector<AggregateId> all_aggregates;
  for (const auto& aggregate_and_rest : demands_) {
    all_aggregates.emplace_back(aggregate_and_rest.first);
  }

  std::shuffle(all_aggregates.begin(), all_aggregates.end(), *rnd);

  std::map<AggregateId, DemandAndFlowCount> new_demands;
  for (size_t i = 0; i < all_aggregates.size(); ++i) {
    const AggregateId& aggregate = all_aggregates[i];
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(demands_, aggregate);

    if (i >= aggregate_count) {
      new_demands[aggregate] = demand_and_flow_count;
      continue;
    }

    nc::net::Bandwidth demand = demand_and_flow_count.first;
    double flow_count = demand_and_flow_count.second;

    nc::net::Bandwidth new_demand = nc::net::Bandwidth::FromBitsPerSecond(
        Pick(demand.bps(), demand_fraction, rnd));

    double new_flow_count = Pick(flow_count, flow_count_fraction, rnd);
    new_demands[aggregate] = {new_demand, new_flow_count};
  }

  auto new_tm = nc::make_unique<TrafficMatrix>(graph_);
  new_tm->demands_ = std::move(new_demands);
  return new_tm;
}

void RoutingConfiguration::AddRouteAndFraction(
    const AggregateId& aggregate_id,
    const std::vector<RouteAndFraction>& routes_and_fractions) {
  CHECK(!nc::ContainsKey(configuration_, aggregate_id));

  // Just to be on the safe side will check that the sum of all fractions is 1
  double total = 0.0;
  for (const auto& route_and_fraction : routes_and_fractions) {
    CHECK(route_and_fraction.second != 0);
    total += route_and_fraction.second;
  }
  CHECK(total <= 1.001 && total >= 0.999) << "Bad total " << total;
  configuration_[aggregate_id] = routes_and_fractions;
}

const std::vector<RouteAndFraction>& RoutingConfiguration::FindRoutesOrDie(
    const AggregateId& aggregate_id) const {
  return nc::FindOrDieNoPrint(configuration_, aggregate_id);
}

std::string RoutingConfiguration::ToString() const {
  std::string base_to_string = TrafficMatrix::ToString();
  std::vector<std::string> out;
  double total = 0;
  for (const auto& aggregate_and_routes : configuration_) {
    const AggregateId& aggregate = aggregate_and_routes.first;
    double aggregate_contribution;
    out.emplace_back(AggregateToString(aggregate, &aggregate_contribution));
    total += aggregate_contribution;
  }

  return nc::StrCat(base_to_string, "\n", nc::Join(out, "\n"), "\ntotal: ",
                    total, "\n");
}

std::string RoutingConfiguration::AggregateToString(
    const AggregateId& aggregate, double* aggregate_contribution) const {
  const std::vector<RouteAndFraction>& routes =
      nc::FindOrDieNoPrint(configuration_, aggregate);
  const DemandAndFlowCount& demand_and_flows =
      nc::FindOrDieNoPrint(demands(), aggregate);

  std::vector<std::string> routes_str;
  double total_contribution = 0;
  for (const RouteAndFraction& route_and_fraction : routes) {
    const nc::net::Walk* path = route_and_fraction.first;
    double fraction = route_and_fraction.second;

    double delay_ms =
        std::chrono::duration<double, std::milli>(path->delay()).count();
    double flow_count = demand_and_flows.second * fraction;
    double contribution = delay_ms * flow_count;
    total_contribution += contribution;

    std::string route_str =
        nc::StrCat(path->ToStringNoPorts(*graph_), " : ", fraction,
                   " (contribution ", contribution, "ms)");
    routes_str.emplace_back(route_str);
  }

  if (aggregate_contribution != nullptr) {
    *aggregate_contribution = total_contribution;
  }

  return nc::StrCat(aggregate.ToString(*graph_), " -> ",
                    nc::Join(routes_str, ","), " (total ", total_contribution,
                    "ms)");
}

void RoutingConfiguration::ToHTML(nc::viz::HtmlPage* out) const {
  using namespace std::chrono;

  OverSubModel model(*this);
  const nc::net::GraphLinkMap<double>& link_to_load = model.link_to_load();
  const std::map<const nc::net::Walk*, nc::net::Bandwidth> path_to_bw =
      model.per_flow_bandwidth_map();

  // Relates each link with the aggregates and paths that cross it. For each
  // path also records the number of flows.
  nc::net::GraphLinkMap<std::map<
      AggregateId, std::vector<std::pair<const nc::net::Walk*, size_t>>>>
      per_link_state;

  std::vector<nc::viz::PathData> paths;
  for (const auto& aggregate_and_routes : configuration_) {
    const AggregateId& aggregate_id = aggregate_and_routes.first;
    const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(demands(), aggregate_id);

    for (const RouteAndFraction& route_and_fraction : routes) {
      std::string label = std::to_string(route_and_fraction.second);
      paths.emplace_back(route_and_fraction.first, label);

      double fraction = route_and_fraction.second;
      const nc::net::Walk* path = route_and_fraction.first;
      for (nc::net::GraphLinkIndex link : path->links()) {
        std::map<AggregateId,
                 std::vector<std::pair<const nc::net::Walk*, size_t>>>&
            aggregates_over_link = per_link_state[link];

        std::vector<std::pair<const nc::net::Walk*, size_t>>& paths_over_link =
            aggregates_over_link[aggregate_id];
        paths_over_link.emplace_back(path,
                                     fraction * demand_and_flow_count.second);
      }
    }
  }

  std::vector<nc::viz::EdgeData> edges;
  for (const auto& link_and_load : link_to_load) {
    nc::net::GraphLinkIndex link = link_and_load.first;
    const nc::net::GraphLink* link_ptr = graph_->GetLink(link);
    double delay_ms =
        std::chrono::duration<double, std::milli>(link_ptr->delay()).count();
    double load = *(link_and_load.second);

    // The tooltip will show the paths that cross the link, grouped by
    // aggregate, and for each path the number of flows and total volume.
    const std::map<AggregateId,
                   std::vector<std::pair<const nc::net::Walk*, size_t>>>&
        aggregates_over_link = per_link_state.GetValueOrDie(link);
    std::string tooltip;
    nc::StrAppend(&tooltip, "Link load: ", load, " out of ",
                  nc::StrCat(link_ptr->bandwidth().Mbps(), "Mbps, delay: ",
                             delay_ms, "ms<br>"));
    for (const auto& aggregate_id_and_paths : aggregates_over_link) {
      const AggregateId& aggregate_id = aggregate_id_and_paths.first;
      nc::StrAppend(&tooltip, aggregate_id.ToString(*graph_), ":<br>");

      for (const auto& path_and_flows : aggregate_id_and_paths.second) {
        const nc::net::Walk* path = path_and_flows.first;
        uint32_t flow_count = path_and_flows.second;
        nc::net::Bandwidth per_flow_bw = nc::FindOrDie(path_to_bw, path);
        nc::net::Bandwidth total_bw = per_flow_bw * flow_count;

        nc::StrAppend(&tooltip, "&nbsp;&nbsp;&nbsp;&nbsp;",
                      path->ToStringIdsOnly(*graph_), " ", flow_count,
                      nc::StrCat(" flows ", total_bw.Mbps(), "Mbps<br>"));
      }
    }

    std::vector<double> loads = {load};
    edges.emplace_back(link, loads, tooltip, delay_ms);
  }

  nc::viz::DisplayMode display_mode("default");
  nc::viz::GraphToHTML(edges, paths, {display_mode}, *graph_, out);
}

std::string TrafficMatrix::ToString() const {
  std::vector<std::string> out;
  for (const auto& aggregate_and_demand : demands_) {
    const AggregateId& aggregate = aggregate_and_demand.first;
    out.emplace_back(AggregateToString(aggregate));
  }

  return nc::Join(out, "\n");
}

std::string TrafficMatrix::AggregateToString(
    const AggregateId& aggregate) const {
  const DemandAndFlowCount& demand_and_flows_count =
      nc::FindOrDieNoPrint(demands_, aggregate);
  nc::net::Bandwidth demand = demand_and_flows_count.first;
  size_t flow_count = demand_and_flows_count.second;

  return nc::StrCat(aggregate.ToString(*graph_), " -> ", demand.Mbps(), "Mbps ",
                    std::to_string(flow_count), " flows");
}

std::string TrafficMatrix::SummaryToString() const {
  std::string graph_summary = graph_->Stats().ToString();
  size_t aggregate_count = demands_.size();

  auto demand_matrix = ToDemandMatrix();
  double mcsf = demand_matrix->MaxCommodityScaleFractor();

  std::vector<nc::net::Bandwidth> demands;
  std::vector<size_t> flow_counts;
  for (const auto& aggregate_and_demand : demands_) {
    const DemandAndFlowCount& demand_and_flow_count =
        aggregate_and_demand.second;
    nc::net::Bandwidth demand = demand_and_flow_count.first;
    size_t flow_count = demand_and_flow_count.second;

    demands.emplace_back(demand);
    flow_counts.emplace_back(flow_count);
  }

  std::vector<nc::net::Bandwidth> demands_percentiles =
      nc::Percentiles(&demands);
  std::string demands_str = nc::Substitute(
      "[min: $0 Mbps, 50p: $1 Mbps, 90p: $2 Mbps, max: $3 Mbps]",
      demands_percentiles[0].Mbps(), demands_percentiles[50].Mbps(),
      demands_percentiles[90].Mbps(), demands_percentiles[100].Mbps());

  std::vector<size_t> flows_percentiles = nc::Percentiles(&flow_counts);
  std::string flow_counts_str = nc::Substitute(
      "[min: $0, 50p: $1, 90p: $2, max: $3]", flows_percentiles[0],
      flows_percentiles[50], flows_percentiles[90], flows_percentiles[100]);

  double fraction_of_all = static_cast<double>(aggregate_count) /
                           (graph_->NodeCount() * (graph_->NodeCount() - 1));
  return nc::Substitute(
      "$0\nnumber of aggregates: $1 ($2%), scale factor: $3\naggregates: "
      "$4\nflow counts: $5\n",
      graph_summary, aggregate_count, fraction_of_all * 100, mcsf, demands_str,
      flow_counts_str);
}

const nc::net::Walk* PickPath(const std::vector<RouteAndFraction>& routes,
                              double init_p) {
  double p = init_p;
  for (const auto& route : routes) {
    double f = route.second;
    if (p < f) {
      return route.first;
    }

    p -= f;
  }

  LOG(FATAL) << "Should not happen";
  return nullptr;
}

static constexpr size_t kTryCount = 1000;

static std::pair<double, double> GetFractionDelta(
    const std::vector<RouteAndFraction>& prev,
    const std::vector<RouteAndFraction>& next) {
  std::mt19937 rnd(1);
  std::uniform_real_distribution<double> dist(0, 1.0);

  double change_count = 0;
  double on_longer_path_count = 0;
  for (size_t i = 0; i < kTryCount; ++i) {
    double p = dist(rnd);
    const nc::net::Walk* p1 = PickPath(prev, p);
    const nc::net::Walk* p2 = PickPath(next, p);

    if (*p1 != *p2) {
      ++change_count;
      if (p2->delay() > p1->delay()) {
        ++on_longer_path_count;
      }
    }
  }

  return {change_count / kTryCount, on_longer_path_count / kTryCount};
}

static void GetRouteCounts(const std::vector<RouteAndFraction>& prev,
                           const std::vector<RouteAndFraction>& next,
                           size_t* add_count, size_t* update_count,
                           size_t* remove_count) {
  for (const RouteAndFraction& prev_fraction : prev) {
    const nc::net::Walk* prev_path = prev_fraction.first;

    bool diff = false;
    bool found = false;
    for (const RouteAndFraction& next_fraction : next) {
      const nc::net::Walk* next_path = next_fraction.first;

      if (*prev_path == *next_path) {
        found = true;

        diff = std::abs(prev_fraction.second - next_fraction.second) > 0.001;
        break;
      }
    }

    if (found) {
      if (diff) {
        ++(*update_count);
      }
    } else {
      ++(*remove_count);
    }
  }

  for (const RouteAndFraction& next_fraction : next) {
    const nc::net::Walk* next_path = next_fraction.first;
    bool found = false;
    for (const RouteAndFraction& prev_fraction : prev) {
      const nc::net::Walk* prev_path = prev_fraction.first;
      if (*prev_path == *next_path) {
        found = true;
      }
    }

    if (!found) {
      ++(*add_count);
    }
  }
}

nc::net::Delay RoutingConfiguration::TotalPerFlowDelay() const {
  double total = 0;
  for (const auto& aggregate_id_and_routes : configuration_) {
    const AggregateId& aggregate_id = aggregate_id_and_routes.first;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(demands(), aggregate_id);

    for (const RouteAndFraction& route_and_fraction :
         aggregate_id_and_routes.second) {
      const nc::net::Walk* path = route_and_fraction.first;
      nc::net::Delay path_delay = path->delay();
      double fraction = route_and_fraction.second;
      double num_flows = demand_and_flow_count.second * fraction;

      total += path_delay.count() * num_flows;
    }
  }

  return nc::net::Delay(static_cast<size_t>(total));
}

std::unique_ptr<RoutingConfiguration> RoutingConfiguration::Copy() const {
  auto out = nc::make_unique<RoutingConfiguration>(
      *static_cast<const TrafficMatrix*>(this));
  out->configuration_ = configuration_;
  return out;
}

RoutingConfigurationDelta RoutingConfiguration::GetDifference(
    const RoutingConfiguration& other) const {
  RoutingConfigurationDelta out;

  size_t total_flow_count = 0;
  double total_flow_count_delta = 0;
  double total_flow_count_delta_longer_path = 0;
  nc::net::Bandwidth total_volume = nc::net::Bandwidth::Zero();
  nc::net::Bandwidth total_volume_delta = nc::net::Bandwidth::Zero();
  nc::net::Bandwidth total_volume_delta_longer_path =
      nc::net::Bandwidth::Zero();
  CHECK(other.configuration_.size() == configuration_.size());
  for (const auto& aggregate_id_and_routes : configuration_) {
    const AggregateId& aggregate_id = aggregate_id_and_routes.first;
    const std::vector<RouteAndFraction>& route_and_fractions_this =
        aggregate_id_and_routes.second;
    const std::vector<RouteAndFraction>& route_and_fractions_other =
        nc::FindOrDieNoPrint(other.configuration_, aggregate_id);

    AggregateDelta& aggregate_delta = out.aggregates[aggregate_id];
    std::tie(aggregate_delta.fraction_delta,
             aggregate_delta.fraction_on_longer_path) =
        GetFractionDelta(route_and_fractions_this, route_and_fractions_other);

    GetRouteCounts(route_and_fractions_this, route_and_fractions_other,
                   &aggregate_delta.routes_added,
                   &aggregate_delta.routes_updated,
                   &aggregate_delta.routes_removed);

    const DemandAndFlowCount& demand_and_flow_count_other =
        nc::FindOrDieNoPrint(other.demands(), aggregate_id);

    size_t flow_count = demand_and_flow_count_other.second;
    nc::net::Bandwidth volume = demand_and_flow_count_other.first;

    total_flow_count_delta += aggregate_delta.fraction_delta * flow_count;
    total_flow_count_delta_longer_path +=
        aggregate_delta.fraction_on_longer_path * flow_count;
    total_flow_count += flow_count;
    total_volume_delta_longer_path +=
        volume * aggregate_delta.fraction_on_longer_path;
    total_volume_delta += volume * aggregate_delta.fraction_delta;
    total_volume += volume;
  }

  nc::net::Delay total_delay = TotalPerFlowDelay();
  nc::net::Delay total_delay_other = TotalPerFlowDelay();
  out.total_per_flow_delay_delta_absolute = total_delay_other - total_delay;
  out.total_per_flow_delay_delta = (total_delay_other - total_delay).count() /
                                   static_cast<double>(total_delay.count());

  out.total_flow_fraction_delta = total_flow_count_delta / total_flow_count;
  out.total_volume_fraction_delta = total_volume_delta / total_volume;
  out.total_flow_fraction_on_longer_path =
      total_flow_count_delta_longer_path / total_flow_count;
  out.total_volume_fraction_on_longer_path =
      total_volume_delta_longer_path / total_volume;
  return out;
}

std::tuple<size_t, size_t, size_t> RoutingConfigurationDelta::TotalRoutes()
    const {
  size_t total_added = 0;
  size_t total_removed = 0;
  size_t total_updated = 0;
  for (const auto& aggregate_and_delta : aggregates) {
    const AggregateDelta& delta = aggregate_and_delta.second;
    total_added += delta.routes_added;
    total_removed += delta.routes_removed;
    total_updated += delta.routes_updated;
  }

  return std::make_tuple(total_added, total_removed, total_updated);
}

nc::net::Bandwidth AggregateHistory::mean_rate() const {
  using namespace std::chrono;
  milliseconds total_time = bins_.size() * bin_size_;
  double total_bytes = 0;
  for (uint64_t bin : bins_) {
    total_bytes += bin;
  }

  double total_time_sec = duration_cast<duration<double>>(total_time).count();
  return nc::net::Bandwidth::FromBitsPerSecond((total_bytes * 8.0) /
                                               total_time_sec);
}

nc::net::Bandwidth AggregateHistory::max_rate() const {
  using namespace std::chrono;

  uint64_t max_bin = *std::max_element(bins_.begin(), bins_.end());
  double total_time_sec = duration_cast<duration<double>>(bin_size_).count();
  return nc::net::Bandwidth::FromBitsPerSecond((max_bin * 8.0) /
                                               total_time_sec);
}

std::chrono::milliseconds AggregateHistory::MaxQueueAtRate(
    nc::net::Bandwidth bandwidth) const {
  double bins_in_second = 1000.0 / bin_size_.count();
  double bytes_per_bin = (bandwidth.bps() / 8.0) / bins_in_second;

  double max_time_to_drain = 0;
  double leftover = 0;
  for (uint64_t bin : bins_) {
    leftover += (bin - bytes_per_bin);
    leftover = std::max(0.0, leftover);

    double time_to_drain_sec = leftover * 8.0 / bandwidth.bps();
    if (max_time_to_drain < time_to_drain_sec) {
      max_time_to_drain = time_to_drain_sec;
    }
  }

  return std::chrono::milliseconds(
      static_cast<uint64_t>(max_time_to_drain * 1000));
}

std::vector<nc::net::Bandwidth> AggregateHistory::PerSecondMeans() const {
  CHECK(1000 % bin_size_.count() == 0);
  size_t bins_in_second = 1000 / bin_size_.count();

  std::vector<nc::net::Bandwidth> out;
  uint64_t current_bytes = 0;
  for (size_t i = 0; i < bins_.size(); ++i) {
    current_bytes += bins_[i];
    if ((i + 1) % bins_in_second == 0) {
      out.emplace_back(
          nc::net::Bandwidth::FromBitsPerSecond(current_bytes * 8));
      current_bytes = 0;
    }
  }

  return out;
}

AggregateHistory AggregateHistory::AddRate(nc::net::Bandwidth rate) const {
  double bins_in_second = 1000.0 / bin_size_.count();
  double bytes_per_bin = (rate.bps() / 8.0) / bins_in_second;
  std::vector<uint64_t> bins_copy = bins_;
  for (uint64_t& bin : bins_copy) {
    bin += bytes_per_bin;
  }

  return {bins_copy, bin_size_, flow_count_};
}

}  // namespace ctr
