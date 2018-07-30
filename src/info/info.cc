#include "info.h"

#include <gflags/gflags.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/map_util.h>
#include <ncode/stats.h>
#include <ncode/strutil.h>
#include <sys/errno.h>
#include <sys/fcntl.h>
#include <unistd.h>
#include <chrono>
#include <cstring>
#include <memory>
#include <string>
#include <utility>

#include "../common.h"
#include "../opt/opt.h"

DEFINE_double(sp_fraction, 1.4,
              "How far from the SP a path can be, for LLPD computation");
DEFINE_double(link_fraction_limit, 0.7,
              "At least this much of the SP's links can be routed around, for "
              "LLPD computation");

namespace ctr {

using namespace std::chrono;

static constexpr size_t kDiscreteMultiplier = 10000;

static double GetLLPDForTopology(const nc::net::GraphStorage& graph) {
  return ctr::GetFractionOfPairsAboveLinkFraction(graph, FLAGS_sp_fraction,
                                                  FLAGS_link_fraction_limit);
}

static info::DiscreteDistribution ToDiscreteDistributionProto(
    const nc::DiscreteDistribution<uint64_t>& distribution,
    double value_multiplier) {
  info::DiscreteDistribution out;
  out.set_value_multiplier(value_multiplier);

  for (const auto& value_and_count : distribution.counts()) {
    uint64_t value = value_and_count.first;
    uint64_t count = value_and_count.second;

    info::ValueAndCount* vc = out.add_values_and_counts();
    vc->set_value(value);
    vc->set_count(count);
  }
  return out;
}

static info::RealDistribution ToRealDistributionProto(
    std::vector<double>* values) {
  info::RealDistribution out;

  std::vector<double> percentiles = nc::Percentiles(values);
  nc::SummaryStats summary_stats;
  for (double value : *values) {
    summary_stats.Add(value);
  }

  out.set_count(summary_stats.count());
  out.set_sum(summary_stats.sum());
  out.set_sum_squared(summary_stats.sum_squared());
  for (double p : percentiles) {
    out.add_percentiles(p);
  }

  return out;
}

static info::RealDistribution ToRealDistributionProto(
    const nc::DiscreteDistribution<uint64_t>& distribution) {
  info::RealDistribution out;

  std::vector<uint64_t> percentiles = distribution.Percentiles();
  nc::SummaryStats summary_stats;
  for (const auto& value_and_count : distribution.counts()) {
    double real_value =
        value_and_count.first / static_cast<double>(kDiscreteMultiplier);
    uint64_t count = value_and_count.second;
    summary_stats.AddCount(real_value, count);
  }

  out.set_count(summary_stats.count());
  out.set_sum(summary_stats.sum());
  out.set_sum_squared(summary_stats.sum_squared());
  for (double p : percentiles) {
    out.add_percentiles(p / static_cast<double>(kDiscreteMultiplier));
  }

  return out;
}

static info::TopologyInfo GenerateTopologyInfo(
    uint64_t id, const nc::net::GraphStorage& graph, const std::string& name) {
  info::TopologyInfo out;
  out.set_topology_id(id);
  out.set_name(name);

  out.set_link_count(graph.LinkCount());
  nc::net::GraphStats stats = graph.Stats();
  out.set_node_count(stats.nodes_count);
  out.set_link_count(stats.links_count);
  out.set_unidirectional_link_count(stats.unidirectional_links);
  out.set_multiple_link_count(stats.multiple_links);
  out.set_diameter_hops(stats.sp_hops.Max());
  out.set_diameter_micros(stats.sp_delays_micros.Max());
  out.set_llpd(GetLLPDForTopology(graph));

  *(out.mutable_node_in_degrees()) =
      ToDiscreteDistributionProto(stats.node_in_degrees, 1.0);
  *(out.mutable_node_out_degrees()) =
      ToDiscreteDistributionProto(stats.node_out_degrees, 1.0);
  *(out.mutable_link_capacities()) =
      ToDiscreteDistributionProto(stats.link_capacities_bps, 1.0);
  *(out.mutable_link_delays()) =
      ToDiscreteDistributionProto(stats.link_delays_micros, 1.0);
  *(out.mutable_shortest_path_delays()) =
      ToDiscreteDistributionProto(stats.sp_delays_micros, 1.0);
  *(out.mutable_shortest_path_hop_counts()) =
      ToDiscreteDistributionProto(stats.sp_hops, 1.0);

  for (nc::net::GraphLinkIndex link : graph.AllLinks()) {
    const nc::net::GraphLink* link_ptr = graph.GetLink(link);
    info::Link* link_info = out.add_links();
    link_info->set_delay_micros(
        duration_cast<microseconds>(link_ptr->delay()).count());
    link_info->set_rate_bps(link_ptr->bandwidth().bps());
    link_info->set_source_index(link_ptr->src());
    link_info->set_destination_index(link_ptr->src());
  }

  for (nc::net::GraphNodeIndex node : graph.AllNodes()) {
    out.add_node_names(graph.GetNode(node)->id());
  }

  return out;
}

static info::TrafficMatrixInfo GenerateTrafficMatrixInfo(
    uint64_t id, uint64_t topology_id, const ctr::TrafficMatrix& tm) {
  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix = tm.ToDemandMatrix();
  info::TrafficMatrixInfo out;
  out.set_traffic_matrix_id(id);
  out.set_topology_id(topology_id);

  uint64_t demand_count = tm.demands().size();
  double node_count = tm.graph()->NodeCount();
  out.set_demand_count(demand_count);
  out.set_demand_fraction(demand_count / (node_count * (node_count - 1)));
  out.set_total_demand_bps(demand_matrix->TotalLoad().bps());
  out.set_is_trivial(demand_matrix->IsTriviallySatisfiable());

  std::string locality_string =
      nc::FindOrDie(demand_matrix->properties(), "locality");
  double locality;
  CHECK(nc::safe_strtod(locality_string, &locality))
      << "Unable to parse locality";
  out.set_locality(locality);
  out.set_max_commodity_scale_factor(
      demand_matrix->MaxCommodityScaleFactor({}, 1.0));

  nc::DiscreteDistribution<uint64_t> demand_rates;
  nc::DiscreteDistribution<uint64_t> demand_flow_counts;
  for (const auto& aggregate_demand_and_flow_count : tm.demands()) {
    const AggregateId& aggregate_id = aggregate_demand_and_flow_count.first;
    const DemandAndFlowCount& demand_and_flow_count =
        aggregate_demand_and_flow_count.second;
    demand_rates.Add(demand_and_flow_count.first.bps());
    demand_flow_counts.Add(demand_and_flow_count.second);

    info::Demand* demand_info = out.add_demands();
    demand_info->set_egress_index(aggregate_id.dst());
    demand_info->set_ingress_index(aggregate_id.src());
    demand_info->set_flow_count(demand_and_flow_count.second);
    demand_info->set_rate_bps(demand_and_flow_count.first.bps());

    nc::net::Delay sp_delay = aggregate_id.GetSPDelay(*tm.graph());
    demand_info->set_shortest_path_delay_micros(
        duration_cast<microseconds>(sp_delay).count());
  }
  *(out.mutable_demand_rates()) =
      ToDiscreteDistributionProto(demand_rates, 1.0);
  *(out.mutable_demand_flow_counts()) =
      ToDiscreteDistributionProto(demand_flow_counts, 1.0);

  return out;
}

static info::RoutingInfo GenerateRoutingInfo(
    const std::string& routing_system, uint64_t id, uint64_t topology_id,
    uint64_t traffic_matrix_id, const ctr::RoutingConfiguration& rc) {
  info::RoutingInfo out;
  out.set_routing_id(id);
  out.set_topology_id(topology_id);
  out.set_traffic_matrix_id(traffic_matrix_id);
  out.set_routing_system(routing_system);

  double total_count = rc.demands().size();
  double overloaded_fraction = rc.OverloadedAggregates().size() / total_count;
  out.set_congested_demands_fraction(overloaded_fraction);

  const nc::net::GraphLinkMap<double>& link_utilizations =
      rc.LinkUtilizations();
  std::vector<double> utilizations;
  nc::net::GraphLinkSet overloaded_links;
  for (const auto& link_and_utilization : link_utilizations) {
    nc::net::GraphLinkIndex link = link_and_utilization.first;
    double utilization = *(link_and_utilization.second);
    utilizations.emplace_back(utilization);
    if (utilization > 1.001) {
      overloaded_links.insert(link);
    }
  }

  *(out.mutable_link_utilizations()) = ToRealDistributionProto(&utilizations);

  nc::DiscreteDistribution<uint64_t> latency_stretch_dist;
  nc::DiscreteDistribution<uint64_t> flow_delay_dist;
  nc::DiscreteDistribution<uint64_t> uncongested_latency_stretch_dist;
  nc::DiscreteDistribution<uint64_t> uncongested_flow_delay_dist;
  double total_delay = 0;
  double total_sp_delay = 0;
  double total_flow_count = 0;
  double total_congested_flow_count = 0;
  for (const auto& aggregate_and_routes : rc.routes()) {
    const AggregateId& id = aggregate_and_routes.first;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(rc.demands(), id);
    size_t aggregate_flow_count = demand_and_flow_count.second;
    total_flow_count += aggregate_flow_count;

    nc::net::Delay sp_delay = id.GetSPDelay(*(rc.graph()));
    uint64_t sp_delay_ms = duration_cast<milliseconds>(sp_delay).count();
    sp_delay_ms = std::max(static_cast<uint64_t>(1), sp_delay_ms);
    total_sp_delay += sp_delay_ms * demand_and_flow_count.second;

    for (const RouteAndFraction& route_and_fraction :
         aggregate_and_routes.second) {
      const nc::net::Walk* path = route_and_fraction.first;
      nc::net::Delay delay = path->delay();
      uint64_t delay_ms = duration_cast<milliseconds>(delay).count();
      delay_ms = std::max(static_cast<uint64_t>(1), delay_ms);

      double fraction = route_and_fraction.second;
      double num_flows = demand_and_flow_count.second * fraction;

      uint64_t num_flows_int = static_cast<uint64_t>(num_flows);
      num_flows_int = std::max(static_cast<uint64_t>(1), num_flows_int);

      double rel_stretch = delay_ms / static_cast<double>(sp_delay_ms);
      bool overloaded = path->ContainsAny(overloaded_links);

      latency_stretch_dist.Add(rel_stretch * kDiscreteMultiplier,
                               num_flows_int);
      flow_delay_dist.Add(delay_ms, num_flows_int);
      if (overloaded) {
        total_congested_flow_count += num_flows;
      } else {
        uncongested_latency_stretch_dist.Add(rel_stretch * kDiscreteMultiplier,
                                             num_flows_int);
        uncongested_flow_delay_dist.Add(delay_ms, num_flows_int);
      }

      total_delay += delay_ms * num_flows;
    }
  }

  double latency_stretch = total_delay / total_sp_delay;
  out.set_total_latency_stretch(latency_stretch);
  out.set_congested_flows_fraction(total_congested_flow_count /
                                   total_flow_count);
  *(out.mutable_flow_delays()) =
      ToDiscreteDistributionProto(flow_delay_dist, 1.0);
  *(out.mutable_flow_delays_uncongested()) =
      ToDiscreteDistributionProto(uncongested_flow_delay_dist, 1.0);
  *(out.mutable_latency_stretch()) =
      ToRealDistributionProto(latency_stretch_dist);
  *(out.mutable_latency_stretch_uncongested()) =
      ToRealDistributionProto(uncongested_latency_stretch_dist);

  return out;
}

void InfoStorage::PrivateAddTopology(const info::TopologyInfo& topology_info) {
  uint64_t id = topology_info.topology_id();
  topology_infos_[id] = topology_info;
}

void InfoStorage::PrivateAddTrafficMatrix(
    const info::TrafficMatrixInfo& tm_info) {
  uint64_t id = tm_info.traffic_matrix_id();
  uint64_t topology_id = tm_info.topology_id();

  info::TrafficMatrixInfo& new_info = traffic_matrix_infos_[id];
  new_info = tm_info;
  topology_to_tm_infos_[topology_id].emplace_back(&new_info);
}

void InfoStorage::PrivateAddRouting(const info::RoutingInfo& routing_info) {
  uint64_t id = routing_info.routing_id();
  uint64_t tm_id = routing_info.traffic_matrix_id();
  const std::string& routing_system = routing_info.routing_system();

  info::RoutingInfo& new_info = routing_infos_[id];
  new_info = routing_info;
  tm_to_routing_infos_[tm_id].emplace_back(&new_info);
  routing_system_to_tm_to_routing_infos_[routing_system][tm_id] = &new_info;
}

uint64_t InfoStorage::AddTopology(const nc::net::GraphStorage& graph,
                                  const std::string& name) {
  std::lock_guard<std::mutex> lock_guard(mu_);
  uint64_t id = GetId();
  PrivateAddTopology(GenerateTopologyInfo(id, graph, name));
  return id;
}

uint64_t InfoStorage::AddTrafficMatrix(uint64_t topology_id,
                                       const ctr::TrafficMatrix& tm) {
  std::lock_guard<std::mutex> lock_guard(mu_);
  uint64_t id = GetId();
  PrivateAddTrafficMatrix(GenerateTrafficMatrixInfo(id, topology_id, tm));
  return id;
}

uint64_t InfoStorage::AddRouting(const std::string& routing_system,
                                 uint64_t topology_id,
                                 uint64_t traffic_matrix_id,
                                 const ctr::RoutingConfiguration& routing) {
  std::lock_guard<std::mutex> lock_guard(mu_);
  uint64_t id = GetId();
  PrivateAddRouting(GenerateRoutingInfo(routing_system, id, topology_id,
                                        traffic_matrix_id, routing));
  return id;
}

const info::TopologyInfo* InfoStorage::GetTopologyInfoOrNull(
    uint64_t topology_id) const {
  return nc::FindOrNull(topology_infos_, topology_id);
}

std::vector<const info::TopologyInfo*> InfoStorage::GetAllTopologyInfos()
    const {
  std::vector<const info::TopologyInfo*> out;
  for (const auto& id_and_info : topology_infos_) {
    out.emplace_back(&(id_and_info.second));
  }
  return out;
}

std::vector<const info::TrafficMatrixInfo*>
InfoStorage::GetAllTrafficMatrixInfos(uint64_t topology_id) const {
  const auto* infos = nc::FindOrNull(topology_to_tm_infos_, topology_id);
  if (infos == nullptr) {
    return {};
  }
  return *infos;
}

const info::TrafficMatrixInfo* InfoStorage::GetTrafficMatrixInfoOrNull(
    uint64_t traffic_matrix_id) const {
  return nc::FindOrNull(traffic_matrix_infos_, traffic_matrix_id);
}

const info::RoutingInfo* InfoStorage::GetRoutingInfoOrNull(
    uint64_t routing_id) const {
  return nc::FindOrNull(routing_infos_, routing_id);
}

std::vector<const info::RoutingInfo*> InfoStorage::GetAllRoutingInfos(
    uint64_t traffic_matrix_id) const {
  const auto* infos = nc::FindOrNull(tm_to_routing_infos_, traffic_matrix_id);
  if (infos == nullptr) {
    return {};
  }
  return *infos;
}

const info::RoutingInfo* InfoStorage::GetRoutingInfoOrNull(
    uint64_t traffic_matrix_id, const std::string& routing_system) const {
  const auto* tm_to_routing_infos =
      nc::FindOrNull(routing_system_to_tm_to_routing_infos_, routing_system);
  if (tm_to_routing_infos == nullptr) {
    return {};
  }

  const auto* info = nc::FindOrNull(*tm_to_routing_infos, traffic_matrix_id);
  if (info == nullptr) {
    return {};
  }

  return *info;
}

uint64_t InfoStorage::GetId() {
  while (true) {
    uint64_t next_id = distribution_(rnd_);
    if (nc::ContainsKey(topology_infos_, next_id) ||
        nc::ContainsKey(traffic_matrix_infos_, next_id) ||
        nc::ContainsKey(routing_infos_, next_id)) {
      continue;
    }

    return next_id;
  }
}

template <typename T>
static bool WriteDelimitedTo(
    const T& entry, google::protobuf::io::FileOutputStream* output_stream) {
  ::google::protobuf::io::CodedOutputStream coded_output(output_stream);
  const int size = entry.ByteSize();
  coded_output.WriteVarint32(size);

  uint8_t* buffer = coded_output.GetDirectBufferForNBytesAndAdvance(size);
  if (buffer != nullptr) {
    // Optimization:  The message fits in one buffer, so use the faster
    // direct-to-array serialization path.
    entry.SerializeWithCachedSizesToArray(buffer);
  } else {
    // Slightly-slower path when the message is multiple buffers.
    entry.SerializeWithCachedSizes(&coded_output);
    if (coded_output.HadError()) {
      return false;
    }
  }

  return true;
}

void InfoStorage::DumpToFile(const std::string& file) const {
  CHECK(!file.empty()) << "Need output file";
  int fd = open(file.c_str(), O_WRONLY | O_TRUNC | O_CREAT,
                S_IREAD | S_IWRITE | S_IRGRP | S_IROTH | S_ISUID);

  CHECK(fd != -1) << "Bad output file " << file;
  google::protobuf::io::FileOutputStream output_stream(fd);

  for (const auto& id_and_info : topology_infos_) {
    info::Info info;
    *(info.mutable_topology_info()) = id_and_info.second;
    CHECK(WriteDelimitedTo(info, &output_stream));
  }
  for (const auto& id_and_info : traffic_matrix_infos_) {
    info::Info info;
    *(info.mutable_traffic_matrix_info()) = id_and_info.second;
    CHECK(WriteDelimitedTo(info, &output_stream));
  }
  for (const auto& id_and_info : routing_infos_) {
    info::Info info;
    *(info.mutable_routing_info()) = id_and_info.second;
    CHECK(WriteDelimitedTo(info, &output_stream));
  }
  output_stream.Close();
}

template <typename T>
static bool ReadDelimitedFrom(T* message,
                              google::protobuf::io::CodedInputStream* input) {
  // Read the size.
  uint32_t size;
  if (!input->ReadVarint32(&size)) {
    return false;
  }

  // Tell the stream not to read beyond that size.
  google::protobuf::io::CodedInputStream::Limit limit = input->PushLimit(size);

  // Parse the message.
  if (!message->MergeFromCodedStream(input)) {
    return false;
  }

  if (!input->ConsumedEntireMessage()) {
    return false;
  }

  // Release the limit.
  input->PopLimit(limit);

  return true;
}

#if defined(__APPLE__) && defined(__MACH__)
#define lseek64 lseek
#define open64 open
#endif

void InfoStorage::ReadFromFile(const std::string& file) {
  int fd = open(file.c_str(), O_RDONLY);
  CHECK(fd > 0) << "Bad input file " << file << ": " << strerror(errno);
  CHECK(lseek64(fd, 0, SEEK_SET) != -1);

  uint64_t topology_count = 0;
  uint64_t tm_count = 0;
  uint64_t routing_count = 0;

  google::protobuf::io::FileInputStream file_input_stream(fd);
  info::Info info;
  while (true) {
    google::protobuf::io::CodedInputStream coded_input(&file_input_stream);
    if (!ReadDelimitedFrom(&info, &coded_input)) {
      break;
    }

    if (info.has_topology_info()) {
      PrivateAddTopology(info.topology_info());
      ++topology_count;
    } else if (info.has_traffic_matrix_info()) {
      PrivateAddTrafficMatrix(info.traffic_matrix_info());
      ++tm_count;
    } else if (info.has_routing_info()) {
      PrivateAddRouting(info.routing_info());
      ++routing_count;
    }

    info.Clear();
  }

  LOG(INFO) << nc::Substitute(
      "Loaded $0 topologies, $1 traffic matrices and $2 routing solutions",
      topology_count, tm_count, routing_count);
}

std::unique_ptr<nc::net::GraphStorage> InfoToGraph(
    const info::TopologyInfo& info) {
  nc::net::GraphBuilder graph_builder;
  for (const auto& link : info.links()) {
    nc::net::Bandwidth bw =
        nc::net::Bandwidth::FromBitsPerSecond(link.rate_bps());
    nc::net::Delay delay = std::chrono::duration_cast<nc::net::Delay>(
        std::chrono::microseconds(link.delay_micros()));
    const std::string& src = info.node_names(link.source_index());
    const std::string& dst = info.node_names(link.destination_index());
    graph_builder.AddLink({src, dst, bw, delay});
  }

  return nc::make_unique<nc::net::GraphStorage>(graph_builder);
}

std::unique_ptr<ctr::TrafficMatrix> InfoToTM(
    const nc::net::GraphStorage& graph, const info::TrafficMatrixInfo& info) {
  auto out = nc::make_unique<ctr::TrafficMatrix>(&graph);
  for (const auto& demand : info.demands()) {
    auto ingress_index = nc::net::GraphNodeIndex(demand.ingress_index());
    auto egress_index = nc::net::GraphNodeIndex(demand.egress_index());

    nc::net::Bandwidth bw =
        nc::net::Bandwidth::FromBitsPerSecond(demand.rate_bps());
    uint64_t flow_count = demand.flow_count();

    out->AddDemand({ingress_index, egress_index}, {bw, flow_count});
  }

  return out;
}

}  // namespace ctr
