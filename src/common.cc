#include "ncode_common/src/common.h"
#include "common.h"

#include <tuple>

#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/strutil.h"

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

void TrafficMatrix::AddDemand(const AggregateId& aggregate_id,
                              const DemandAndFlowCount& demand_and_flow_count) {
  CHECK(!nc::ContainsKey(demands_, aggregate_id));
  demands_[aggregate_id] = demand_and_flow_count;
}

TrafficMatrix::TrafficMatrix(
    const nc::lp::DemandMatrix& demand_matrix,
    const std::map<nc::lp::SrcAndDst, double>& flow_counts)
    : graph_(demand_matrix.graph()) {
  for (const auto& element : demand_matrix.elements()) {
    double flow_count =
        nc::FindWithDefault(flow_counts, {element.src, element.dst}, 1.0);
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
    double demand_fraction, double flow_count_fraction,
    std::mt19937* rnd) const {
  std::map<AggregateId, DemandAndFlowCount> new_demands;
  for (const auto& aggregate_and_demand : demands_) {
    const AggregateId& aggregate = aggregate_and_demand.first;
    nc::net::Bandwidth demand = aggregate_and_demand.second.first;
    double flow_count = aggregate_and_demand.second.second;

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
  CHECK(total <= 1.001 && total >= 0.999);
  configuration_[aggregate_id] = routes_and_fractions;
}

const std::vector<RouteAndFraction>& RoutingConfiguration::FindRoutesOrDie(
    const AggregateId& aggregate_id) const {
  return nc::FindOrDieNoPrint(configuration_, aggregate_id);
}

std::string RoutingConfiguration::ToString(
    const nc::net::GraphStorage& graph) const {
  std::vector<std::string> out;
  for (const auto& aggregate_and_routes : configuration_) {
    const AggregateId& aggregate = aggregate_and_routes.first;
    const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;

    std::vector<std::string> routes_str;
    for (const RouteAndFraction& route_and_fraction : routes) {
      std::string route_str =
          nc::StrCat(route_and_fraction.first->ToStringNoPorts(graph), ":",
                     route_and_fraction.second);
      routes_str.emplace_back(route_str);
    }

    out.emplace_back(
        nc::StrCat(aggregate.ToString(graph), "->", nc::Join(routes_str, ",")));
  }

  return nc::Join(out, "\n");
}

}  // namespace ctr
