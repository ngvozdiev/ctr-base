#include "info_server.h"

#include <gflags/gflags.h>
#include <google/protobuf/repeated_field.h>
#include <ncode/common.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/net/net_common.h>
#include <ncode/substitute.h>
#include <ncode/viz/server.h>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#include "../common.h"
#include "tm_gen.h"

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

#define EVALUATE_NUM_FIELD(name)                    \
  if (mask.name()) {                                \
    if ((info.name() < constraint.lower_limit()) || \
        (info.name() > constraint.upper_limit())) { \
      return false;                                 \
    }                                               \
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
  SET_COMPLEX_FIELD(node_names);
  SET_SIMPLE_FIELD(name);
  SET_SIMPLE_FIELD(date);
  SET_SIMPLE_FIELD(description);
  SET_DISCRETE_DIST(shortest_path_delays);
  SET_DISCRETE_DIST(shortest_path_hop_counts);

  return out;
}

static bool ConstraintMatches(const info::TopologyInfo& info,
                              const info::FieldConstraint& constraint) {
  if (constraint.has_topology_mask()) {
    const info::TopologyInfoMask& mask = constraint.topology_mask();
    EVALUATE_NUM_FIELD(node_count);
    EVALUATE_NUM_FIELD(link_count);
    EVALUATE_NUM_FIELD(unidirectional_link_count);
    EVALUATE_NUM_FIELD(multiple_link_count);
    EVALUATE_NUM_FIELD(diameter_hops);
    EVALUATE_NUM_FIELD(diameter_micros);
    EVALUATE_NUM_FIELD(llpd);
  }

  return true;
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
  SET_SIMPLE_FIELD(max_commodity_scale_factor);
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

static bool ConstraintMatches(const info::TrafficMatrixInfo& info,
                              const info::FieldConstraint& constraint) {
  if (constraint.has_traffic_matrix_mask()) {
    const info::TrafficMatrixInfoMask& mask = constraint.traffic_matrix_mask();
    EVALUATE_NUM_FIELD(demand_count);
    EVALUATE_NUM_FIELD(demand_fraction);
    EVALUATE_NUM_FIELD(total_demand_bps);
    EVALUATE_NUM_FIELD(locality);
    EVALUATE_NUM_FIELD(max_commodity_scale_factor);
  }

  return true;
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

static bool ConstraintMatches(const info::RoutingInfo& info,
                              const info::FieldConstraint& constraint) {
  if (constraint.has_routing_mask()) {
    const info::RoutingInfoMask& mask = constraint.routing_mask();
    EVALUATE_NUM_FIELD(congested_flows_fraction);
    EVALUATE_NUM_FIELD(congested_demands_fraction);
    EVALUATE_NUM_FIELD(total_latency_stretch);
  }

  return true;
}

template <typename T>
static bool ConstraintMatchesExpression(
    const T& element,
    const info::ConstraintOrExpression& constraint_or_expression) {
  if (constraint_or_expression.has_constraint()) {
    return ConstraintMatches(element, constraint_or_expression.constraint());
  }

  if (constraint_or_expression.has_expression()) {
    const info::ConstraintExpression& expression =
        constraint_or_expression.expression();
    switch (expression.type()) {
      case info::ConstraintExpression::NOT:
        CHECK(expression.has_op_one());
        return !ConstraintMatchesExpression(element, expression.op_one());
      case info::ConstraintExpression::AND:
        CHECK(expression.has_op_one());
        CHECK(expression.has_op_two());
        return ConstraintMatchesExpression(element, expression.op_one()) &&
               ConstraintMatchesExpression(element, expression.op_two());
      case info::ConstraintExpression::OR:
        CHECK(expression.has_op_one());
        CHECK(expression.has_op_two());
        return ConstraintMatchesExpression(element, expression.op_one()) ||
               ConstraintMatchesExpression(element, expression.op_two());
      default:
        LOG(ERROR) << "Invalid expression type " << expression.type();
        return false;
    }
  }

  return true;
}

template <typename T>
static void Filter(std::vector<T>* values, std::function<bool(T)> filter) {
  std::vector<T> out;
  for (const auto& v : *values) {
    if (filter(v)) {
      out.emplace_back(v);
    }
  }

  *values = out;
}

template <typename T>
static void Filter(std::set<T>* values, std::function<bool(T)> filter) {
  std::set<T> out;
  for (const auto& v : *values) {
    if (filter(v)) {
      out.insert(v);
    }
  }

  *values = out;
}

InfoServer::InfoServer(std::unique_ptr<InfoStorage> info_storage)
    : nc::ProtobufServer<ctr::info::Request, ctr::info::Response>(
          nc::viz::TCPServerConfig(FLAGS_info_server_port, 1 << 10),
          FLAGS_info_server_threads),
      info_storage_(std::move(info_storage)) {}

std::unique_ptr<ctr::info::Response> InfoServer::HandleSelect(
    const info::SelectInfoRequest& request) const {
  std::vector<const info::TopologyInfo*> selected_topologies;
  std::set<const info::TrafficMatrixInfo*> selected_traffic_matrices;
  std::set<const info::RoutingInfo*> selected_routings;

  bool return_topologies = false;
  bool return_tms = false;
  bool return_routing = false;
  if (request.topology_mask().SerializeAsString() !=
      info::TopologyInfoMask::default_instance().SerializeAsString()) {
    return_topologies = true;
  }
  if (request.traffic_matrix_mask().SerializeAsString() !=
      info::TrafficMatrixInfoMask::default_instance().SerializeAsString()) {
    return_topologies = true;
    return_tms = true;
  }
  if (request.routing_mask().SerializeAsString() !=
      info::RoutingInfoMask::default_instance().SerializeAsString()) {
    return_topologies = true;
    return_tms = true;
    return_routing = true;
  }

  if (return_topologies) {
    if (request.topology_id_size() == 0) {
      for (const info::TopologyInfo* info :
           info_storage_->GetAllTopologyInfos()) {
        selected_topologies.emplace_back(info);
      }
    } else {
      for (uint64_t topology_id : request.topology_id()) {
        const info::TopologyInfo* info =
            info_storage_->GetTopologyInfoOrNull(topology_id);
        if (info != nullptr) {
          selected_topologies.emplace_back(info);
        }
      }
    }

    Filter<const info::TopologyInfo*>(
        &selected_topologies, [&request](const info::TopologyInfo* info) {
          return ConstraintMatchesExpression(*info, request.constraints());
        });
  }

  if (return_tms) {
    if (request.traffic_matrix_id_size() == 0) {
      for (const info::TopologyInfo* topology : selected_topologies) {
        std::vector<const info::TrafficMatrixInfo*> infos =
            info_storage_->GetAllTrafficMatrixInfos(topology->topology_id());
        selected_traffic_matrices.insert(infos.begin(), infos.end());
      }
    } else {
      for (uint64_t traffic_matrix_id : request.traffic_matrix_id()) {
        const info::TrafficMatrixInfo* info =
            info_storage_->GetTrafficMatrixInfoOrNull(traffic_matrix_id);
        if (info != nullptr) {
          selected_traffic_matrices.insert(info);
        }
      }
    }

    Filter<const info::TrafficMatrixInfo*>(
        &selected_traffic_matrices,
        [&request](const info::TrafficMatrixInfo* info) {
          return ConstraintMatchesExpression(*info, request.constraints());
        });
  }

  if (return_routing) {
    if (request.routing_id_size() == 0 && request.routing_systems_size() == 0) {
      for (const info::TrafficMatrixInfo* tm_info : selected_traffic_matrices) {
        std::vector<const info::RoutingInfo*> infos =
            info_storage_->GetAllRoutingInfos(tm_info->traffic_matrix_id());
        selected_routings.insert(infos.begin(), infos.end());
      }
    } else {
      for (uint64_t routing_id : request.routing_id()) {
        const info::RoutingInfo* info =
            info_storage_->GetRoutingInfoOrNull(routing_id);
        if (info != nullptr) {
          selected_routings.insert(info);
        }
      }

      for (const std::string& routing_system : request.routing_systems()) {
        for (const info::TrafficMatrixInfo* tm_info :
             selected_traffic_matrices) {
          const info::RoutingInfo* info = info_storage_->GetRoutingInfoOrNull(
              tm_info->traffic_matrix_id(), routing_system);
          if (info != nullptr) {
            selected_routings.insert(info);
          }
        }
      }
    }

    Filter<const info::RoutingInfo*>(
        &selected_routings, [&request](const info::RoutingInfo* info) {
          return ConstraintMatchesExpression(*info, request.constraints());
        });
  }

  auto out = nc::make_unique<ctr::info::Response>();
  info::SelectInfoResponse* select_response =
      out->mutable_select_info_response();
  for (const info::TopologyInfo* topology_info : selected_topologies) {
    *(select_response->add_topology_info()) =
        ApplyTopologyMask(*topology_info, request.topology_mask());
  }
  for (const info::TrafficMatrixInfo* traffic_matrix_info :
       selected_traffic_matrices) {
    *(select_response->add_traffic_matrix_info()) = ApplyTrafficInfoMask(
        *traffic_matrix_info, request.traffic_matrix_mask());
  }
  for (const info::RoutingInfo* routing_info : selected_routings) {
    *(select_response->add_routing_info()) =
        ApplyRoutingInfoMask(*routing_info, request.routing_mask());
  }

  return out;
}

std::unique_ptr<ctr::info::Response> InfoServer::HandleRequest(
    const ctr::info::Request& request) {
  LOG(INFO) << "R " << request.DebugString();
  if (request.has_select_info_request()) {
    return HandleSelect(request.select_info_request());
  }
  if (request.has_traffic_matrix_generate_request()) {
    return HandleTMGenRequest(request.traffic_matrix_generate_request());
  }

  LOG(INFO) << "Empty request";
  return nc::make_unique<ctr::info::Response>();
}

static std::unique_ptr<ctr::info::Response> GetTMGenResponseWithError(
    const std::string& error) {
  auto response = nc::make_unique<ctr::info::Response>();
  response->mutable_generate_traffic_matrix_response()->set_error(error);
  return response;
}

std::unique_ptr<ctr::info::Response> InfoServer::HandleTMGenRequest(
    const info::TrafficMatrixGenerateRequest& request) const {
  if (request.flow_count_distribution_case() ==
      info::TrafficMatrixGenerateRequest::FLOW_COUNT_DISTRIBUTION_NOT_SET) {
    return GetTMGenResponseWithError(
        "Need one of constant_flow_count or demand_based_total_flow_count");
  }

  uint64_t topology_id = request.topology_id();
  const info::TopologyInfo* topology_info =
      info_storage_->GetTopologyInfoOrNull(topology_id);
  if (topology_info == nullptr) {
    return GetTMGenResponseWithError(
        nc::Substitute("Unable to find topology with ID $0", topology_id));
  }

  std::unique_ptr<nc::net::GraphStorage> graph = InfoToGraph(*topology_info);
  std::mt19937 gen(request.seed());
  auto demand_matrix =
      GenerateRoughan(*graph, nc::net::Bandwidth::FromGBitsPerSecond(1), &gen);

  demand_matrix = LocalizeDemandMatrix(*demand_matrix, request.locality());
  demand_matrix = BalanceReverseDemands(*demand_matrix, 0.1);
  double csf = demand_matrix->MaxCommodityScaleFactor({}, 1.0);
  CHECK(csf != 0);
  CHECK(csf == csf);
  demand_matrix = demand_matrix->Scale(csf);

  double load = 1.0 / request.max_commodity_scale_factor();
  demand_matrix = demand_matrix->Scale(load);

  std::unique_ptr<ctr::TrafficMatrix> tm;
  switch (request.flow_count_distribution_case()) {
    case info::TrafficMatrixGenerateRequest::kConstantFlowCount: {
      tm = TrafficMatrix::ConstantFromDemandMatrix(
          *demand_matrix, request.constant_flow_count());
      break;
    }

    case info::TrafficMatrixGenerateRequest::kDemandBasedTotalFlowCount: {
      tm = TrafficMatrix::DistributeFromDemandMatrix(
          *demand_matrix, request.demand_based_total_flow_count());
      break;
    }

    case info::TrafficMatrixGenerateRequest::FLOW_COUNT_DISTRIBUTION_NOT_SET: {
      LOG(FATAL) << "Should not happen";
    }
  }

  uint64_t traffic_matrix_id =
      info_storage_->AddTrafficMatrix(topology_id, *tm);
  auto response = nc::make_unique<ctr::info::Response>();
  response->mutable_generate_traffic_matrix_response()->set_traffic_matrix_id(
      traffic_matrix_id);
  return response;
}

void InfoServer::LogMessage(const LoggedMessage& logged_message) {
  info::RequestLogEntry log_entry;
  *(log_entry.mutable_request()) = *(logged_message.request);
  log_entry.set_request_size_bytes(logged_message.request->ByteSize());
  log_entry.set_response_size_bytes(logged_message.response_size_bytes);

  const nc::viz::TCPConnectionInfo& connection_info =
      logged_message.header_and_message->tcp_connection_info;
  log_entry.set_connection_id(connection_info.connection_id);
  log_entry.set_remote_ip(connection_info.remote_ip);
  log_entry.set_remote_port(connection_info.remote_port);

  log_entry.set_request_rx_ms(
      logged_message.header_and_message->time_rx.count());
  log_entry.set_request_started_processing_ms(
      logged_message.processing_started.count());
  log_entry.set_request_done_processing_ms(
      logged_message.processing_done.count());
  log_entry.set_response_sent_ms(logged_message.response_sent.count());
  log_entry.set_worker_index(logged_message.worker_index);

  LOG(INFO) << log_entry.DebugString();
}

}  // namespace ctr
