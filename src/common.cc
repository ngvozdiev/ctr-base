#include "common.h"

#include <ncode/common.h>
#include <ncode/event_queue.h>
#include <ncode/free_list.h>
#include <ncode/htsim/packet.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/lp/mc_flow.h>
#include <ncode/net/algorithm.h>
#include <ncode/net/net_common.h>
#include <ncode/perfect_hash.h>
#include <ncode/stats.h>
#include <ncode/strutil.h>
#include <ncode/substitute.h>
#include <ncode/viz/graph.h>
#include <stddef.h>
#include <initializer_list>
#include <random>
#include <ratio>
#include <set>
#include <string>
#include <tuple>

#include "opt/oversubscription_model.h"
#include "opt/path_provider.h"

namespace ctr {

constexpr char RoutingConfiguration::kDefaultOptimizerString[];

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
  return GetSP(graph)->delay();
}

std::unique_ptr<nc::net::Walk> AggregateId::GetSP(
    const nc::net::GraphStorage& graph) const {
  std::unique_ptr<nc::net::Walk> sp =
      nc::net::ShortestPathWithConstraints(src_, dst_, graph, {});
  CHECK(sp);
  return sp;
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

std::unique_ptr<TrafficMatrix> TrafficMatrix::ProportionalFromDemandMatrix(
    const nc::lp::DemandMatrix& demand_matrix,
    size_t top_aggregate_flow_count) {
  nc::net::Bandwidth max_demand = nc::net::Bandwidth::Zero();
  for (const auto& element : demand_matrix.elements()) {
    max_demand = std::max(max_demand, element.demand);
  }

  // Will assign flow counts so that the max bandwidth aggregate has
  // kTopAggregateFlowCount, and all other aggregates proportionally less.
  std::map<AggregateId, DemandAndFlowCount> demands_and_counts;
  for (const auto& element : demand_matrix.elements()) {
    size_t flow_count =
        top_aggregate_flow_count * (element.demand / max_demand);
    flow_count = std::max(1ul, flow_count);

    AggregateId id(element.src, element.dst);
    demands_and_counts[id] = {element.demand, flow_count};
  }

  return nc::make_unique<TrafficMatrix>(demand_matrix.graph(),
                                        demands_and_counts);
}

std::unique_ptr<TrafficMatrix> TrafficMatrix::DistributeFromDemandMatrix(
    const nc::lp::DemandMatrix& demand_matrix, size_t total_flow_count) {
  nc::net::Bandwidth total_demand = nc::net::Bandwidth::Zero();
  for (const auto& element : demand_matrix.elements()) {
    total_demand += element.demand;
  }

  std::map<AggregateId, DemandAndFlowCount> demands_and_counts;
  for (const auto& element : demand_matrix.elements()) {
    double f = element.demand / total_demand;
    size_t flow_count = total_flow_count * f;
    flow_count = std::max(1ul, flow_count);

    AggregateId id(element.src, element.dst);
    demands_and_counts[id] = {element.demand, flow_count};
  }

  return nc::make_unique<TrafficMatrix>(demand_matrix.graph(),
                                        demands_and_counts);
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

std::pair<nc::net::Bandwidth, nc::net::Bandwidth>
TrafficMatrix::MinMaxAggregates() const {
  nc::net::Bandwidth min_bw = nc::net::Bandwidth::Max();
  nc::net::Bandwidth max_bw = nc::net::Bandwidth::Zero();

  for (const auto& aggregate_and_demand : demands_) {
    const DemandAndFlowCount& demand_and_flow_count =
        aggregate_and_demand.second;
    min_bw = std::min(min_bw, demand_and_flow_count.first);
    max_bw = std::max(max_bw, demand_and_flow_count.first);
  }

  return {min_bw, max_bw};
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

std::unique_ptr<TrafficMatrix> TrafficMatrix::AddToDemands(
    nc::net::Bandwidth value, const std::set<AggregateId>& to_extend) const {
  std::map<AggregateId, DemandAndFlowCount> new_demands;
  for (const auto& aggregate_and_demand : demands_) {
    const AggregateId& aggregate = aggregate_and_demand.first;
    const DemandAndFlowCount& demand_and_flow_count =
        aggregate_and_demand.second;
    if (!to_extend.empty() && !nc::ContainsKey(to_extend, aggregate)) {
      new_demands[aggregate] = demand_and_flow_count;
      continue;
    }

    new_demands[aggregate] = {demand_and_flow_count.first + value,
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

void CompetingAggregates::AddAggregatesAndCapacity(
    const AggregateId& aggregate_id,
    const std::vector<AggregatesAndCapacity>& aggregates_and_capacity) {
  CHECK(!nc::ContainsKey(aggregates_, aggregate_id));
  aggregates_[aggregate_id] = aggregates_and_capacity;
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

std::map<AggregateId, std::vector<size_t>> RoutingConfiguration::GetKValues()
    const {
  std::map<AggregateId, std::vector<size_t>> out;
  size_t i = 0;
  for (const auto& id_and_routes : configuration_) {
    const AggregateId& id = id_and_routes.first;
    std::vector<size_t>& indices = out[id];

    nc::net::KShortestPathsGenerator ksp_gen(id.src(), id.dst(), *graph_, {});
    for (const auto& route : id_and_routes.second) {
      const nc::net::Walk& to_look_for = *route.first;
      indices.emplace_back(ksp_gen.KforPath(to_look_for));
    }

    ++i;
  }

  return out;
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

  for (nc::net::GraphLinkIndex link_index : graph_->AllLinks()) {
    per_link_state[link_index] = {};
  }

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
  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix = ToDemandMatrix();
  std::vector<double> sp_utilizations;
  nc::net::GraphLinkMap<double> per_link_utilization =
      demand_matrix->SPUtilization();
  for (const auto& link_and_utilization : per_link_utilization) {
    sp_utilizations.emplace_back(*link_and_utilization.second);
  }
  std::sort(sp_utilizations.begin(), sp_utilizations.end());

  out.emplace_back(
      nc::StrCat("TM with ", static_cast<uint64_t>(demands_.size()),
                 " demands, scale factor ",
                 demand_matrix->MaxCommodityScaleFactor({}, 1.0),
                 " sp link utilizations: ", nc::Join(sp_utilizations, ",")));
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
  double mcsf = demand_matrix->MaxCommodityScaleFactor({}, 1.0);

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

std::vector<FlowPathChange> GetFractionDeltas(
    const std::vector<RouteAndFraction>& src,
    const std::vector<RouteAndFraction>& dst) {
  // This assumes that each of prev and next are already sorted by delay and
  // their fractions both sum up to 1.
  std::vector<FlowPathChange> out;

  double src_f_cumulative = 0.0;
  double dst_f_cumulative = 0.0;
  double fraction = 0.0;
  size_t src_index = 0;
  size_t dst_index = 0;
  while (true) {
    if (src_index == src.size()) {
      for (size_t i = dst_index; i < dst.size(); ++i) {
        double dst_f = dst_f_cumulative + dst[dst_index].second;
        out.push_back({dst_f - fraction, src.back().first, dst[i].first});
        dst_f_cumulative = dst_f;
        fraction = dst_f;
      }

      break;
    }

    if (dst_index == dst.size()) {
      for (size_t i = src_index; i < src.size(); ++i) {
        double src_f = src_f_cumulative + src[src_index].second;
        out.push_back({src_f - fraction, src[i].first, dst.back().first});
        src_f_cumulative = src_f;
        fraction = src_f;
      }

      break;
    }

    double src_f = src_f_cumulative + src[src_index].second;
    double dst_f = dst_f_cumulative + dst[dst_index].second;
    if (src_f > dst_f) {
      out.push_back(
          {dst_f - fraction, src[src_index].first, dst[dst_index].first});
      ++dst_index;
      dst_f_cumulative = dst_f;
      fraction = dst_f;
    } else if (src_f < dst_f) {
      out.push_back(
          {src_f - fraction, src[src_index].first, dst[dst_index].first});
      ++src_index;
      src_f_cumulative = src_f;
      fraction = src_f;
    } else {
      out.push_back(
          {src_f - fraction, src[src_index].first, dst[dst_index].first});
      ++src_index;
      ++dst_index;
      dst_f_cumulative = dst_f;
      src_f_cumulative = src_f;
      fraction = src_f;
    }
  }

  return out;
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

nc::net::Delay RoutingConfiguration::TotalPerFlowDelay(bool sp) const {
  double total = 0;
  for (const auto& aggregate_id_and_routes : configuration_) {
    const AggregateId& aggregate_id = aggregate_id_and_routes.first;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(demands(), aggregate_id);

    if (sp) {
      nc::net::Delay sp_delay = aggregate_id.GetSPDelay(*graph_);
      total += sp_delay.count() * demand_and_flow_count.second;
      continue;
    }

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

nc::net::GraphLinkMap<double> RoutingConfiguration::LinkUtilizations() const {
  // A map from a link to the total load over the link.
  std::map<nc::net::GraphLinkIndex, nc::net::Bandwidth> link_to_total_load;
  for (const auto& aggregate_and_aggregate_output : configuration_) {
    const AggregateId& aggregate_id = aggregate_and_aggregate_output.first;
    const std::vector<RouteAndFraction>& routes =
        aggregate_and_aggregate_output.second;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(demands(), aggregate_id);
    nc::net::Bandwidth total_aggregate_demand = demand_and_flow_count.first;

    for (const auto& route : routes) {
      const nc::net::Walk* path = route.first;
      double fraction = route.second;
      CHECK(fraction > 0);

      for (nc::net::GraphLinkIndex link : path->links()) {
        link_to_total_load[link] += total_aggregate_demand * fraction;
      }
    }
  }

  nc::net::GraphLinkMap<double> out;
  for (const auto& link_and_total_load : link_to_total_load) {
    nc::net::GraphLinkIndex link_index = link_and_total_load.first;
    const nc::net::GraphLink* link = graph_->GetLink(link_index);

    nc::net::Bandwidth total_load = link_and_total_load.second;
    double utilization = total_load / link->bandwidth();
    out[link_index] = utilization;
  }

  return out;
}

size_t RoutingConfiguration::OverloadedAggregates() const {
  std::map<nc::net::GraphLinkIndex, nc::net::Bandwidth> link_to_total_load;
  nc::net::GraphLinkMap<std::set<AggregateId>> link_to_aggregate;
  for (const auto& aggregate_and_aggregate_output : configuration_) {
    const AggregateId& aggregate_id = aggregate_and_aggregate_output.first;
    const std::vector<RouteAndFraction>& routes =
        aggregate_and_aggregate_output.second;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(demands(), aggregate_id);
    nc::net::Bandwidth total_aggregate_demand = demand_and_flow_count.first;

    for (const auto& route : routes) {
      const nc::net::Walk* path = route.first;
      double fraction = route.second;
      CHECK(fraction > 0);

      for (nc::net::GraphLinkIndex link : path->links()) {
        link_to_total_load[link] += total_aggregate_demand * fraction;
        link_to_aggregate[link].emplace(aggregate_id);
      }
    }
  }

  std::set<AggregateId> overloaded_aggregates;
  for (const auto& link_and_total_load : link_to_total_load) {
    nc::net::GraphLinkIndex link_index = link_and_total_load.first;
    const nc::net::GraphLink* link = graph_->GetLink(link_index);

    nc::net::Bandwidth total_load = link_and_total_load.second;
    double utilization = total_load / link->bandwidth();
    if (utilization > 1) {
      const std::set<AggregateId>& aggregates_crossing_link =
          link_to_aggregate.GetValueOrDie(link_index);
      overloaded_aggregates.insert(aggregates_crossing_link.begin(),
                                   aggregates_crossing_link.end());
    }
  }

  return overloaded_aggregates.size();
}

double RoutingConfiguration::MaxLinkUtilization() const {
  nc::net::GraphLinkMap<double> link_utilizations = LinkUtilizations();
  double max_utilization = 0;
  for (const auto& link_and_utilization : link_utilizations) {
    max_utilization = std::max(max_utilization, *link_and_utilization.second);
  }
  return max_utilization;
}

std::unique_ptr<RoutingConfiguration> RoutingConfiguration::Copy() const {
  auto out = nc::make_unique<RoutingConfiguration>(
      *static_cast<const TrafficMatrix*>(this));
  out->configuration_ = configuration_;
  out->time_to_compute_ = time_to_compute_;
  out->optimizer_string_ = optimizer_string_;
  return out;
}

RoutingConfigurationDelta RoutingConfiguration::GetDifference(
    const RoutingConfiguration& other) const {
  const nc::net::GraphStorage* other_graph = other.graph();
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

    // The other's graph may be different. Need to create a different id.
    const std::string& src_id = graph_->GetNode(aggregate_id.src())->id();
    const std::string& dst_id = graph_->GetNode(aggregate_id.dst())->id();
    AggregateId other_aggregate_id(other_graph->NodeFromStringOrDie(src_id),
                                   other_graph->NodeFromStringOrDie(dst_id));

    const std::vector<RouteAndFraction>& route_and_fractions_this =
        aggregate_id_and_routes.second;
    const std::vector<RouteAndFraction>& route_and_fractions_other =
        nc::FindOrDieNoPrint(other.configuration_, other_aggregate_id);

    AggregateDelta& aggregate_delta = out.aggregates[aggregate_id];
    aggregate_delta.changes =
        GetFractionDeltas(route_and_fractions_this, route_and_fractions_other);

    GetRouteCounts(route_and_fractions_this, route_and_fractions_other,
                   &aggregate_delta.routes_added,
                   &aggregate_delta.routes_updated,
                   &aggregate_delta.routes_removed);

    const DemandAndFlowCount& demand_and_flow_count_other =
        nc::FindOrDieNoPrint(other.demands(), other_aggregate_id);

    size_t flow_count = demand_and_flow_count_other.second;
    nc::net::Bandwidth volume = demand_and_flow_count_other.first;

    total_flow_count_delta += aggregate_delta.FractionDelta() * flow_count;
    total_flow_count_delta_longer_path +=
        aggregate_delta.FractionOnLongerPath() * flow_count;
    total_flow_count += flow_count;
    total_volume_delta_longer_path +=
        volume * aggregate_delta.FractionOnLongerPath();
    total_volume_delta += volume * aggregate_delta.FractionDelta();
    total_volume += volume;
  }

  nc::net::Delay total_delay = TotalPerFlowDelay();
  nc::net::Delay total_delay_other = other.TotalPerFlowDelay();
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

std::string RoutingConfigurationDelta::ToString(
    const nc::net::GraphStorage& graph) const {
  using namespace std::chrono;
  std::string out;
  nc::SubstituteAndAppend(
      &out,
      "flow fraction delta: $0\nvolume fraction delta: $1\nflow fraction on "
      "longer path: $2\nvolume fraction on longer path: $3\ntotal per flow "
      "delay fraction: $4\ntotal per flow delay abs: $5 ms\n",
      total_flow_fraction_delta, total_volume_fraction_delta,
      total_flow_fraction_on_longer_path, total_volume_fraction_on_longer_path,
      total_per_flow_delay_delta,
      duration_cast<milliseconds>(total_per_flow_delay_delta_absolute).count());

  nc::StrAppend(&out, "Aggregates with change:\n");
  for (const auto& aggregate_and_delta : aggregates) {
    const AggregateId& aggregate = aggregate_and_delta.first;
    const AggregateDelta& delta = aggregate_and_delta.second;
    if (delta.FractionDelta() == 0) {
      continue;
    }

    nc::SubstituteAndAppend(
        &out, "$0: fraction delta: $1, fraction on longer path: $2\n",
        aggregate.ToString(graph), delta.FractionDelta(),
        delta.FractionOnLongerPath());
  }

  return out;
}

double AggregateDelta::FractionDelta(double relative_change_threshold) const {
  double total_same = 0;
  for (const auto& path_change : changes) {
    uint64_t from_delay = path_change.from->delay().count();
    uint64_t to_delay = path_change.to->delay().count();

    double relative_change =
        (std::max(from_delay, to_delay) - std::min(from_delay, to_delay)) /
        static_cast<double>(from_delay);
    if (relative_change > relative_change_threshold) {
      continue;
    }

    total_same += path_change.fraction;
  }

  return 1 - total_same;
}

double AggregateDelta::FractionOnLongerPath(
    double relative_change_threshold) const {
  double total_longer = 0;
  for (const auto& path_change : changes) {
    uint64_t from_delay = path_change.from->delay().count();
    uint64_t to_delay = path_change.to->delay().count();
    if (from_delay >= to_delay) {
      continue;
    }

    double relative_change =
        (to_delay - from_delay) / static_cast<double>(from_delay);
    if (relative_change <= relative_change_threshold) {
      continue;
    }

    total_longer += path_change.fraction;
  }

  return total_longer;
}

AggregateHistory::AggregateHistory(nc::net::Bandwidth rate, size_t bin_count,
                                   std::chrono::milliseconds bin_size,
                                   size_t flow_count)
    : bin_size_(bin_size), flow_count_(flow_count) {
  using namespace std::chrono;

  double duration_sec = duration_cast<duration<double>>(bin_size_).count();
  double rate_Bps = rate.bps() / 8.0;
  uint64_t bin_size_bytes = duration_sec * rate_Bps;
  bins_.resize(bin_count, bin_size_bytes);
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

nc::net::Bandwidth AggregateHistory::MaxRateAtQueue(
    std::chrono::milliseconds max_queue) const {
  nc::net::Bandwidth min_bound = mean_rate();
  nc::net::Bandwidth max_bound = max_rate();

  nc::net::Bandwidth curr_estimate = max_bound;
  while (true) {
    CHECK(max_bound >= min_bound);
    nc::net::Bandwidth delta = max_bound - min_bound;
    if (delta <= nc::net::Bandwidth::FromKBitsPerSecond(1)) {
      break;
    }

    nc::net::Bandwidth guess = min_bound + (max_bound - min_bound) / 2;
    std::chrono::milliseconds queue_at_guess = MaxQueueAtRate(guess);
    bool is_ok = queue_at_guess < max_queue;
    if (is_ok) {
      curr_estimate = guess;
      // Will try to go lower.
      max_bound = guess;
    } else {
      min_bound = guess;
    }
  }

  return curr_estimate;
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

AggregateHistory AggregateHistory::SubtractRate(nc::net::Bandwidth rate) const {
  double bins_in_second = 1000.0 / bin_size_.count();
  double bytes_per_bin = (rate.bps() / 8.0) / bins_in_second;
  std::vector<uint64_t> bins_copy = bins_;
  for (uint64_t& bin : bins_copy) {
    bin -= bytes_per_bin;
  }

  return {bins_copy, bin_size_, flow_count_};
}

TLDRRequest::TLDRRequest(
    nc::net::IPAddress ip_src, nc::net::IPAddress ip_dst,
    nc::EventQueueTime time_sent, uint64_t request_id, bool quick,
    const std::map<AggregateId, AggregateHistory>& aggregates)
    : nc::htsim::Message(ip_src, ip_dst, kTLDRRequestType, time_sent),
      round_id_(request_id),
      aggregates_(aggregates),
      quick_(quick) {}

const std::map<AggregateId, AggregateHistory>& TLDRRequest::aggregates() {
  return aggregates_;
}

std::string TLDRRequest::ToString() const {
  std::string out =
      nc::Substitute("MSG $0 -> $1 : TLDRRequest id $2 ",
                     nc::net::IPToStringOrDie(five_tuple_.ip_src()),
                     nc::net::IPToStringOrDie(five_tuple_.ip_dst()), round_id_);
  for (const auto& id_and_rate : aggregates_) {
    AggregateHistory rate_and_flow_count = id_and_rate.second;
    uint32_t bin_count = rate_and_flow_count.bins().size();
    nc::StrAppend(&out, nc::StrCat("{", bin_count, "bins, q:", quick_, "}"));
  }

  nc::StrAppend(&out, " sent at ", time_sent().Raw());
  return out;
}

nc::htsim::PacketPtr TLDRRequest::Duplicate() const {
  return nc::GetFreeList<TLDRRequest>().New(five_tuple_.ip_src(),
                                            five_tuple_.ip_dst(), time_sent(),
                                            round_id_, quick_, aggregates_);
}

TLDRForceRequest::TLDRForceRequest(nc::net::IPAddress ip_src,
                                   nc::net::IPAddress ip_dst,
                                   nc::EventQueueTime time_sent)
    : nc::htsim::Message(ip_src, ip_dst, kTLDRForceRequestType, time_sent) {}

std::string TLDRForceRequest::ToString() const {
  return nc::Substitute("MSG $0 -> $1 : TLDRForceRequest",
                        nc::net::IPToStringOrDie(five_tuple_.ip_src()),
                        nc::net::IPToStringOrDie(five_tuple_.ip_dst()));
}

nc::htsim::PacketPtr TLDRForceRequest::Duplicate() const {
  return nc::GetFreeList<TLDRForceRequest>().New(
      five_tuple_.ip_src(), five_tuple_.ip_dst(), time_sent());
}

TLDRTriggerReoptimize::TLDRTriggerReoptimize(nc::net::IPAddress ip_src,
                                             nc::net::IPAddress ip_dst,
                                             nc::EventQueueTime time_sent)
    : nc::htsim::Message(ip_src, ip_dst, kTLDRTriggerReoptimizeType,
                         time_sent) {}

std::string TLDRTriggerReoptimize::ToString() const {
  return nc::Substitute("MSG $0 -> $1 : TLDRTriggerReoptimize",
                        nc::net::IPToStringOrDie(five_tuple_.ip_src()),
                        nc::net::IPToStringOrDie(five_tuple_.ip_dst()));
}

nc::htsim::PacketPtr TLDRTriggerReoptimize::Duplicate() const {
  return nc::GetFreeList<TLDRTriggerReoptimize>().New(
      five_tuple_.ip_src(), five_tuple_.ip_dst(), time_sent());
}

TLDRUpdate::TLDRUpdate(
    nc::net::IPAddress ip_src, nc::net::IPAddress ip_dst,
    nc::EventQueueTime time_sent,
    const std::map<AggregateId, AggregateUpdateState>& aggregate_updates,
    const std::map<AggregateId, AggregateHistory>& additional_histories)
    : ::nc::htsim::Message(ip_src, ip_dst, kTLDRUpdateType, time_sent),
      aggregates_(aggregate_updates),
      additional_histories_(additional_histories) {}

nc::htsim::PacketPtr TLDRUpdate::Duplicate() const {
  return nc::GetFreeList<TLDRUpdate>().New(five_tuple_.ip_src(),
                                           five_tuple_.ip_dst(), time_sent(),
                                           aggregates_, additional_histories_);
}

std::string TLDRUpdate::ToString() const {
  std::string out =
      nc::Substitute("MSG $0 -> $1 : TLDRUpdate",
                     nc::net::IPToStringOrDie(five_tuple_.ip_src()),
                     nc::net::IPToStringOrDie(five_tuple_.ip_dst()));
  return out;
}

ctr::AggregateHistory GetDummyHistory(nc::net::Bandwidth rate,
                                      std::chrono::milliseconds bin_size,
                                      std::chrono::milliseconds duration,
                                      size_t flow_count) {
  size_t bin_size_ms = bin_size.count();
  size_t init_window_ms = duration.count();
  size_t bins_count = init_window_ms / bin_size_ms;

  double bins_in_second = 1000.0 / bin_size_ms;
  CHECK(bins_in_second > 0);
  double bytes_per_bin = (rate.bps() / 8.0) / bins_in_second;
  CHECK(bytes_per_bin > 0);

  std::vector<uint64_t> bins;
  for (size_t i = 0; i < bins_count; ++i) {
    bins.emplace_back(bytes_per_bin);
  }

  return {bins, std::chrono::milliseconds(bin_size_ms), flow_count};
}

static std::vector<nc::net::GraphNodeIndex> NodesFromPath(
    const nc::net::Walk& path, const nc::net::GraphStorage& graph) {
  std::vector<nc::net::GraphNodeIndex> out;
  for (nc::net::GraphLinkIndex link : path.links()) {
    nc::net::GraphNodeIndex node = graph.GetLink(link)->src();
    out.emplace_back(node);
  }

  nc::net::GraphLinkIndex last_link = path.links().back();
  nc::net::GraphNodeIndex last_node = graph.GetLink(last_link)->dst();
  out.emplace_back(last_node);
  return out;
}

std::string RoutingConfiguration::SerializeToText(
    const std::vector<std::string>& node_order) const {
  std::map<std::string, uint32_t> indices;
  for (size_t i = 0; i < node_order.size(); ++i) {
    indices[node_order[i]] = i;
  }

  std::string out;
  nc::StrAppend(&out, optimizer_string_, ",", time_to_compute_.count(), "\n");
  for (const auto& aggregate_and_paths : configuration_) {
    for (const auto& path_and_fraction : aggregate_and_paths.second) {
      const nc::net::Walk* path = path_and_fraction.first;
      double fraction = path_and_fraction.second;

      std::vector<nc::net::GraphNodeIndex> nodes =
          NodesFromPath(*path, *graph_);

      std::vector<std::string> strings;
      for (nc::net::GraphNodeIndex node : nodes) {
        const std::string& node_id = graph_->GetNode(node)->id();
        strings.emplace_back(nc::StrCat(nc::FindOrDie(indices, node_id)));
      }

      strings.emplace_back(nc::StrCat(fraction));
      nc::StrAppend(&out, nc::Join(strings, ","), "\n");
    }
  }

  return out;
}

std::unique_ptr<RoutingConfiguration>
RoutingConfiguration::LoadFromSerializedText(
    const TrafficMatrix& base_matrix,
    const std::vector<std::string>& node_order, const std::string& text,
    PathProvider* path_provider) {
  std::vector<std::string> lines = nc::Split(text, "\n");
  CHECK(!lines.empty());

  std::vector<std::string> first_line = nc::Split(lines[0], ",");
  CHECK(first_line.size() == 2);

  std::string optimizer_string = first_line[0];
  uint32_t time_to_compute;
  CHECK(nc::safe_strtou32(first_line[1], &time_to_compute));

  std::map<AggregateId, std::vector<RouteAndFraction>> routes;
  const nc::net::GraphStorage* graph = base_matrix.graph();

  for (size_t i = 1; i < lines.size(); ++i) {
    const std::string& line = lines[i];
    std::vector<std::string> line_split = nc::Split(line, ",");
    CHECK(line_split.size() >= 3);

    double fraction;
    CHECK(nc::safe_strtod(line_split.back(), &fraction));

    std::vector<nc::net::GraphNodeIndex> nodes_on_path;
    for (size_t j = 0; j < line_split.size() - 1; ++j) {
      uint32_t node_index;
      CHECK(nc::safe_strtou32(line_split[j], &node_index));
      CHECK(node_index < node_order.size());

      const std::string& node_id = node_order[node_index];
      nc::net::GraphNodeIndex node = graph->NodeFromStringOrDie(node_id);
      nodes_on_path.emplace_back(node);
    }
    CHECK(nodes_on_path.size() >= 2);

    nc::net::Links links_on_path;
    for (size_t j = 0; j < nodes_on_path.size() - 1; ++j) {
      auto link = graph->LinkOrDie(nodes_on_path[j], nodes_on_path[j + 1]);
      links_on_path.emplace_back(link);
    }

    auto new_path = nc::make_unique<nc::net::Walk>(links_on_path, *graph);
    const nc::net::Walk* path_ptr =
        path_provider->TakeOwnership(std::move(new_path));

    AggregateId id(nodes_on_path.front(), nodes_on_path.back());
    CHECK(nc::ContainsKey(base_matrix.demands(), id));
    routes[id].emplace_back(path_ptr, fraction);
  }

  auto out = nc::make_unique<RoutingConfiguration>(base_matrix);
  out->time_to_compute_ = std::chrono::milliseconds(time_to_compute);
  out->optimizer_string_ = optimizer_string;

  for (const auto& aggregate_id_and_routes : routes) {
    const AggregateId& id = aggregate_id_and_routes.first;
    out->AddRouteAndFraction(id, aggregate_id_and_routes.second);
  }

  return out;
}

}  // namespace ctr
