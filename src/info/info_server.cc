#include "info_server.h"

#include <gflags/gflags.h>
#include <google/protobuf/repeated_field.h>
#include <ncode/common.h>
#include <ncode/logging.h>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#define SET_SIMPLE_FIELD(name)       \
  if (mask.name()) {                 \
    out.set_##name(original.name()); \
  }

#define SET_COMPLEX_FIELD(name)                \
  if (mask.name()) {                           \
    *(out.mutable_##name()) = original.name(); \
  }

#define SET_DISCRETE_DIST(name)                                               \
  if (mask.name##_multiplier() != 0.0) {                                      \
    *(out.mutable_##name()) =                                                 \
        ScaleDiscreteDistribution(original.name(), mask.name##_multiplier()); \
  }

DEFINE_uint64(info_server_port, 8080, "The port to run the server on");
DEFINE_uint64(info_server_threads, 10,
              "Number of threads to serve requests on");

namespace ctr {

static info::DiscreteDistribution ScaleDiscreteDistribution(
    const info::DiscreteDistribution& original, double multiplier) {
  CHECK(original.value_multiplier() == 1.0);
  std::map<uint64_t, uint64_t> new_value_counts;
  for (const auto& value_and_count : original.values_and_counts()) {
    uint64_t value = value_and_count.value();
    uint64_t count = value_and_count.count();

    uint64_t new_value = static_cast<uint64_t>(value * multiplier);
    new_value_counts[new_value] += count;
  }

  info::DiscreteDistribution out;
  out.set_value_multiplier(multiplier);
  for (const auto& new_value_and_count : new_value_counts) {
    info::ValueAndCount* value_and_count = out.add_values_and_counts();
    value_and_count->set_value(new_value_and_count.first);
    value_and_count->set_count(new_value_and_count.second);
  }

  return out;
}

static info::TopologyInfo ApplyTopologyMask(
    const info::TopologyInfo& original, const info::TopologyInfoMask& mask) {
  info::TopologyInfo out;
  out.set_topology_id(original.topology_id());
  SET_SIMPLE_FIELD(node_count);
  SET_SIMPLE_FIELD(link_count);
  SET_SIMPLE_FIELD(unidirectional_link_count);
  SET_SIMPLE_FIELD(multiple_link_count);
  SET_SIMPLE_FIELD(diameter_hops);
  SET_SIMPLE_FIELD(diameter_micros);
  SET_SIMPLE_FIELD(llpd);
  SET_COMPLEX_FIELD(node_in_degrees);
  SET_COMPLEX_FIELD(node_out_degrees);
  SET_COMPLEX_FIELD(link_capacities);
  SET_COMPLEX_FIELD(link_delays);
  SET_COMPLEX_FIELD(links);
  SET_SIMPLE_FIELD(name);
  SET_SIMPLE_FIELD(date);
  SET_SIMPLE_FIELD(description);
  SET_DISCRETE_DIST(shortest_path_delays);
  SET_DISCRETE_DIST(shortest_path_hop_counts);

  return out;
}

static info::Demand ApplyDemandMask(const info::Demand& original,
                                    const info::DemandMask& mask) {
  info::Demand out;
  SET_SIMPLE_FIELD(ingress_index);
  SET_SIMPLE_FIELD(egress_index);
  SET_SIMPLE_FIELD(rate_bps);
  SET_SIMPLE_FIELD(flow_count);
  SET_SIMPLE_FIELD(shortest_path_delay_micros);
  return out;
}

static info::TrafficMatrixInfo ApplyTrafficInfoMask(
    const info::TrafficMatrixInfo& original,
    const info::TrafficMatrixInfoMask& mask) {
  info::TrafficMatrixInfo out;
  out.set_topology_id(original.topology_id());
  out.set_traffic_matrix_id(original.traffic_matrix_id());
  SET_SIMPLE_FIELD(demand_count);
  SET_SIMPLE_FIELD(demand_fraction);
  SET_SIMPLE_FIELD(total_demand_bps);
  SET_SIMPLE_FIELD(locality);
  SET_SIMPLE_FIELD(limiting_multiplier);
  SET_SIMPLE_FIELD(limiting_rate);
  SET_SIMPLE_FIELD(is_trivial);
  SET_COMPLEX_FIELD(demands_out_degree);
  SET_COMPLEX_FIELD(demands_in_degree);
  SET_DISCRETE_DIST(demand_rates);
  SET_DISCRETE_DIST(demand_flow_counts);

  if (mask.demand_mask().ByteSize() != 0) {
    for (const auto& demand : original.demands()) {
      *(out.add_demands()) = ApplyDemandMask(demand, mask.demand_mask());
    }
  }

  return out;
}

static info::RoutingInfo ApplyRoutingInfoMask(
    const info::RoutingInfo& original, const info::RoutingInfoMask& mask) {
  info::RoutingInfo out;
  out.set_routing_id(original.routing_id());
  out.set_traffic_matrix_id(original.traffic_matrix_id());
  out.set_routing_system(original.routing_system());
  SET_SIMPLE_FIELD(congested_flows_fraction);
  SET_SIMPLE_FIELD(congested_demands_fraction);
  SET_SIMPLE_FIELD(total_latency_stretch);
  SET_COMPLEX_FIELD(link_utilizations);
  SET_COMPLEX_FIELD(latency_stretch);
  SET_COMPLEX_FIELD(latency_stretch_uncongested);
  SET_DISCRETE_DIST(flow_delays);
  SET_DISCRETE_DIST(flow_delays_uncongested);
  return out;
}

InfoServer::InfoServer(std::unique_ptr<InfoStorage> info_storage)
    : nc::ProtobufServer<ctr::info::InfoRequest, ctr::info::InfoResponse>(
          FLAGS_info_server_port, FLAGS_info_server_threads),
      info_storage_(std::move(info_storage)) {}

std::unique_ptr<ctr::info::InfoResponse> InfoServer::HandleRequest(
    const ctr::info::InfoRequest& request) {
  if (request.has_topology_info_request()) {
    return HandleTopologyInfo(request.topology_info_request());
  }
  if (request.has_traffic_info_request()) {
    return HandleTrafficMatrixInfo(request.traffic_info_request());
  }
  if (request.has_routing_info_request()) {
    return HandleRoutingInfo(request.routing_info_request());
  }

  LOG(INFO) << "Empty request";
  return nc::make_unique<ctr::info::InfoResponse>();
}

std::unique_ptr<ctr::info::InfoResponse> InfoServer::HandleTopologyInfo(
    const info::TopologyInfoRequest& request) const {
  std::vector<const info::TopologyInfo*> infos_to_serve;
  if (request.topology_ids_size() == 0) {
    infos_to_serve = info_storage_->GetAllTopologyInfos();
  } else {
    for (uint64_t topology_id : request.topology_ids()) {
      const info::TopologyInfo* info =
          info_storage_->GetTopologyInfoOrNull(topology_id);
      if (info != nullptr) {
        infos_to_serve.emplace_back(info);
      }
    }
  }

  auto out = nc::make_unique<ctr::info::InfoResponse>();
  if (infos_to_serve.empty()) {
    return out;
  }

  info::TopologyInfoResponse* response = out->mutable_topology_info_response();
  for (const info::TopologyInfo* info : infos_to_serve) {
    *(response->add_topology_info()) = ApplyTopologyMask(*info, request.mask());
  }

  return out;
}

std::unique_ptr<ctr::info::InfoResponse> InfoServer::HandleTrafficMatrixInfo(
    const info::TrafficMatrixInfoRequest& request) const {
  std::vector<const info::TrafficMatrixInfo*> infos_to_serve;
  if (request.traffic_matrix_ids_size() == 0) {
    infos_to_serve =
        info_storage_->GetAllTrafficMatrixInfos(request.topology_id());
  } else {
    for (uint64_t topology_id : request.traffic_matrix_ids()) {
      const info::TrafficMatrixInfo* info =
          info_storage_->GetTrafficMatrixInfoOrNull(topology_id);
      if (info != nullptr) {
        infos_to_serve.emplace_back(info);
      }
    }
  }

  auto out = nc::make_unique<ctr::info::InfoResponse>();
  if (infos_to_serve.empty()) {
    return out;
  }

  info::TrafficMatrixInfoResponse* response =
      out->mutable_traffic_info_response();
  for (const info::TrafficMatrixInfo* info : infos_to_serve) {
    *(response->add_traffic_info()) =
        ApplyTrafficInfoMask(*info, request.mask());
  }

  return out;
}

std::unique_ptr<ctr::info::InfoResponse> InfoServer::HandleRoutingInfo(
    const info::RoutingInfoRequest& request) const {
  std::set<const info::RoutingInfo*> infos_to_serve;
  if (request.routing_systems_size() == 0 && request.routing_ids_size()) {
    std::vector<const info::RoutingInfo*> infos_vector =
        info_storage_->GetAllRoutingInfos(request.traffic_matrix_id());
    infos_to_serve.insert(infos_vector.begin(), infos_vector.end());
  } else {
    for (uint64_t routing_id : request.routing_ids()) {
      const info::RoutingInfo* info =
          info_storage_->GetRoutingInfoOrNull(routing_id);
      if (info != nullptr) {
        infos_to_serve.emplace(info);
      }
    }

    for (const std::string& routing_system : request.routing_systems()) {
      const info::RoutingInfo* info = info_storage_->GetRoutingInfoOrNull(
          request.traffic_matrix_id(), routing_system);
      if (info != nullptr) {
        infos_to_serve.emplace(info);
      }
    }
  }

  auto out = nc::make_unique<ctr::info::InfoResponse>();
  if (infos_to_serve.empty()) {
    return out;
  }

  info::RoutingInfoResponse* response = out->mutable_routing_info_response();
  for (const info::RoutingInfo* info : infos_to_serve) {
    *(response->add_routing_info()) =
        ApplyRoutingInfoMask(*info, request.mask());
  }

  return out;
}
}  // namespace ctr
