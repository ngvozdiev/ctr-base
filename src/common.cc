#include "ncode_common/src/common.h"
#include "common.h"

#include <tuple>

#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/net/algorithm.h"

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
  for (const auto& aggregate_and_routes : configuration_) {
    const AggregateId& aggregate = aggregate_and_routes.first;
    const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;

    std::vector<std::string> routes_str;
    for (const RouteAndFraction& route_and_fraction : routes) {
      std::string route_str =
          nc::StrCat(route_and_fraction.first->ToStringNoPorts(*graph_), ":",
                     route_and_fraction.second);
      routes_str.emplace_back(route_str);
    }

    out.emplace_back(nc::StrCat(aggregate.ToString(*graph_), " -> ",
                                nc::Join(routes_str, ",")));
  }

  return nc::StrCat(base_to_string, "\n", nc::Join(out, "\n"));
}

std::string TrafficMatrix::ToString() const {
  std::vector<std::string> out;
  for (const auto& aggregate_and_demand : demands_) {
    const AggregateId& aggregate = aggregate_and_demand.first;
    const DemandAndFlowCount& demand_and_flows_count =
        aggregate_and_demand.second;

    nc::net::Bandwidth demand = demand_and_flows_count.first;
    size_t flow_count = demand_and_flows_count.second;

    out.emplace_back(nc::StrCat(aggregate.ToString(*graph_), " -> ",
                                demand.Mbps(), "Mbps ",
                                std::to_string(flow_count), " flows"));
  }

  return nc::Join(out, "\n");
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

static double GetFractionDelta(const std::vector<RouteAndFraction>& lhs,
                               const std::vector<RouteAndFraction>& rhs) {
  // Will get the difference as 1 - fraction that stays on the same path.
  double total = 0;
  for (const RouteAndFraction& lhs_fraction : lhs) {
    const nc::net::Walk* lhs_path = lhs_fraction.first;
    for (const RouteAndFraction& rhs_fraction : rhs) {
      const nc::net::Walk* rhs_path = rhs_fraction.first;

      if (*lhs_path == *rhs_path) {
        double f_delta = std::min(lhs_fraction.second, rhs_fraction.second);
        total += f_delta;
      }
    }
  }

  return 1 - total;
}

static double GetWeightedPathStretch(
    const std::vector<RouteAndFraction>& routes, nc::net::Delay sp_delay) {
  double sp_delay_sec = std::chrono::duration<double>(sp_delay).count();

  double total = 0;
  for (const RouteAndFraction& route_and_fraction : routes) {
    const nc::net::Walk* path = route_and_fraction.first;
    double path_delay_sec =
        std::chrono::duration<double>(path->delay()).count();
    double stretch = path_delay_sec - sp_delay_sec;
    CHECK(stretch >= 0);

    total += stretch * route_and_fraction.second;
  }

  return total;
}

static double GetPathStretchGain(const std::vector<RouteAndFraction>& lhs,
                                 const std::vector<RouteAndFraction>& rhs,
                                 nc::net::Delay sp_delay) {
  return GetWeightedPathStretch(lhs, sp_delay) -
         GetWeightedPathStretch(rhs, sp_delay);
}

RoutingConfigurationDelta RoutingConfiguration::GetDifference(
    const RoutingConfiguration& other) const {
  RoutingConfigurationDelta out;

  size_t total_flow_count = 0;
  double total_flow_count_delta = 0;
  nc::net::Bandwidth total_volume = nc::net::Bandwidth::Zero();
  nc::net::Bandwidth total_volume_delta = nc::net::Bandwidth::Zero();
  CHECK(other.configuration_.size() == configuration_.size());
  for (const auto& aggregate_id_and_routes : configuration_) {
    const AggregateId& aggregate_id = aggregate_id_and_routes.first;
    const std::vector<RouteAndFraction>& route_and_fractions_this =
        aggregate_id_and_routes.second;
    const std::vector<RouteAndFraction>& route_and_fractions_other =
        nc::FindOrDieNoPrint(other.configuration_, aggregate_id);

    AggregateDelta& aggregate_delta = out.aggregates[aggregate_id];
    aggregate_delta.fraction_delta =
        GetFractionDelta(route_and_fractions_this, route_and_fractions_other);
    aggregate_delta.path_stretch_gain =
        GetPathStretchGain(route_and_fractions_this, route_and_fractions_other,
                           aggregate_id.GetSPDelay(*graph_));

    const DemandAndFlowCount& demand_and_flow_count_this =
        nc::FindOrDieNoPrint(demands(), aggregate_id);
    const DemandAndFlowCount& demand_and_flow_count_other =
        nc::FindOrDieNoPrint(other.demands(), aggregate_id);

    size_t flow_count = std::max(demand_and_flow_count_this.second,
                                 demand_and_flow_count_other.second);
    nc::net::Bandwidth volume = std::max(demand_and_flow_count_this.first,
                                         demand_and_flow_count_other.first);
    total_flow_count_delta += aggregate_delta.fraction_delta * flow_count;
    total_flow_count += flow_count;
    total_volume_delta += volume * aggregate_delta.fraction_delta;
    total_volume += volume;
  }

  out.total_flow_fraction_delta = total_flow_count_delta / total_flow_count;
  out.total_volume_fraction_delta = total_volume_delta / total_volume;
  return out;
}

}  // namespace ctr
