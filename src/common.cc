#include "common.h"

#include <utility>
#include <vector>

#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
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
                              const DemandAndPriority& demand_and_priority) {
  CHECK(!nc::ContainsKey(demands_, aggregate_id));
  demands_[aggregate_id] = demand_and_priority;
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

std::unique_ptr<TrafficMatrix> TrafficMatrix::Randomize(
    double fraction, std::mt19937* rnd) const {
  std::map<AggregateId, DemandAndPriority> new_demands;
  for (const auto& aggregate_and_demand : demands_) {
    const AggregateId& aggregate = aggregate_and_demand.first;
    nc::net::Bandwidth demand = aggregate_and_demand.second.first;

    double demand_mbps = demand.Mbps();
    double delta = demand_mbps * fraction;

    double low = demand_mbps - delta;
    double high = demand_mbps + delta;
    std::uniform_real_distribution<double> dist(low, high);
    double new_demand = dist(*rnd);

    new_demands[aggregate] = {
        nc::net::Bandwidth::FromMBitsPerSecond(new_demand),
        aggregate_and_demand.second.second};
  }

  auto new_tm = nc::make_unique<TrafficMatrix>(graph_);
  new_tm->demands_ = std::move(new_demands);
  return new_tm;
}

// Parses a line of the form <tag> <count> and returns count.
static uint32_t ParseCountOrDie(const std::string& tag,
                                const std::string& line) {
  std::vector<std::string> line_split = nc::Split(line, " ");
  CHECK(line_split.size() == 2);
  CHECK(line_split[0] == tag);

  uint32_t count;
  CHECK(nc::safe_strtou32(line_split[1], &count));
  return count;
}

std::unique_ptr<TrafficMatrix> TrafficMatrix::LoadRepetitaOrDie(
    const std::string& matrix_string,
    const std::vector<std::string>& node_names,
    const nc::net::GraphStorage* graph) {
  std::vector<std::string> all_lines = nc::Split(matrix_string, "\n");
  auto it = all_lines.begin();

  const std::string& demands_line = *it;
  uint32_t num_demands = ParseCountOrDie("DEMANDS", demands_line);

  // Skip free form line.
  ++it;

  std::map<std::pair<nc::net::GraphNodeIndex, nc::net::GraphNodeIndex>, double>
      total_demands;
  for (uint32_t i = 0; i < num_demands; ++i) {
    ++it;

    std::vector<std::string> line_split = nc::Split(*it, " ");
    CHECK(line_split.size() == 4) << *it << " demand " << i;

    uint32_t src_index;
    uint32_t dst_index;
    double demand_kbps;

    CHECK(nc::safe_strtou32(line_split[1], &src_index));
    CHECK(nc::safe_strtou32(line_split[2], &dst_index));
    CHECK(nc::safe_strtod(line_split[3], &demand_kbps));

    CHECK(src_index < node_names.size()) << src_index << " line " << *it;
    CHECK(dst_index < node_names.size()) << dst_index << " line " << *it;

    nc::net::GraphNodeIndex src =
        graph->NodeFromStringOrDie(node_names[src_index]);
    nc::net::GraphNodeIndex dst =
        graph->NodeFromStringOrDie(node_names[dst_index]);

    total_demands[{src, dst}] += demand_kbps;
  }

  std::map<AggregateId, DemandAndPriority> new_demands;
  for (const auto& src_and_dst_and_demand : total_demands) {
    double demand_kbps = src_and_dst_and_demand.second;
    if (demand_kbps < 1) {
      continue;
    }

    nc::net::GraphNodeIndex src;
    nc::net::GraphNodeIndex dst;
    std::tie(src, dst) = src_and_dst_and_demand.first;

    nc::net::Bandwidth demand =
        nc::net::Bandwidth::FromKBitsPerSecond(demand_kbps);
    new_demands[{src, dst}] = demand;
  }

  auto new_tm = nc::make_unique<TrafficMatrix>(graph);
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
