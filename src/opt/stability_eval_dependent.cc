#include "../common.h"

namespace ctr {

std::unique_ptr<TrafficMatrix> ScaleBasedOnOutput(
    const RoutingConfiguration& output,
    const std::map<AggregateId, nc::net::Bandwidth>& capacity_on_sp) {
  for (const auto& aggregate_and_routes : output.routes()) {
    const AggregateId& aggregate_id = aggregate_and_routes.first;

    // Need to figure out what fraction of the aggregate's volume is on
  }
}

}  // namespace ctr
