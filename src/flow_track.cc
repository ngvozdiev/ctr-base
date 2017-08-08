#include "flow_track.h"

#include <stddef.h>
#include <algorithm>
#include <iterator>
#include <string>
#include <unordered_map>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/strutil.h"

namespace nc {

using namespace nc::net;
using namespace nc::pcap;

std::string DataCluster::ToString() const {
  nc::pcap::Timestamp duration = last_byte_at - first_byte_at;
  uint64_t duration_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
  return nc::StrCat("duration: ", duration_ms, "ms, pkts: ", packets);
}

static std::vector<uint64_t> Combine(
    const UnidirectionalTCPFlowState& forward,
    const UnidirectionalTCPFlowState& reverse) {
  std::vector<uint64_t> timestamps_forward = forward.timestamps.Restore();
  std::vector<uint64_t> timestamps_reverse = reverse.timestamps.Restore();

  std::vector<uint64_t> all_values;
  std::move(timestamps_forward.begin(), timestamps_forward.end(),
            std::back_inserter(all_values));
  std::move(timestamps_reverse.begin(), timestamps_reverse.end(),
            std::back_inserter(all_values));

  std::sort(all_values.begin(), all_values.end());
  return all_values;
}

static std::vector<DataCluster> BreakDownStatic(
    Timestamp max_delta, const UnidirectionalTCPFlowState& forward,
    const UnidirectionalTCPFlowState& reverse) {
  std::vector<uint64_t> timestamps_combined = Combine(forward, reverse);

  std::vector<DataCluster> clusters;
  for (size_t i = 0; i < timestamps_combined.size(); ++i) {
    uint64_t time = timestamps_combined[i];
    if (clusters.empty()) {
      clusters.emplace_back(Timestamp(time));
    } else {
      uint64_t prev_time = clusters.back().last_byte_at.count();
      CHECK(time >= prev_time);
      uint64_t delta = time - prev_time;
      if (delta > static_cast<uint64_t>(max_delta.count())) {
        clusters.emplace_back(Timestamp(time));
      }
    }

    DataCluster& current_cluster = clusters.back();
    current_cluster.last_byte_at = Timestamp(time);
    ++current_cluster.packets;
  }

  return clusters;
}

std::vector<DataCluster> BidirectionalTCPFlowState::BreakDown(
    Timestamp max_delta) const {
  return BreakDownStatic(max_delta, *client_to_server_flow_,
                         *server_to_client_flow_);
}

std::string BidirectionalTCPFlowState::ToString() const {
  std::string out = StrCat("Client to server 5-tuple: ",
                           client_to_server_tuple_.ToString(), "\n");
  std::vector<DataCluster> clusters = BreakDown(std::chrono::milliseconds(100));
  std::vector<std::string> cluster_strings;
  for (const DataCluster& cluster : clusters) {
    cluster_strings.emplace_back(cluster.ToString());
  }

  return StrCat(out, Join(cluster_strings, ","), "\n");
}

void FlowTracker::HandleTCP(Timestamp timestamp, const IPHeader& ip_header,
                            const TCPHeader& tcp_header, uint16_t payload_len) {
  nc::Unused(payload_len);

  net::IPAddress src_address(ntohl(ip_header.ip_src.s_addr));
  net::IPAddress dst_address(ntohl(ip_header.ip_dst.s_addr));
  net::AccessLayerPort src_port(ntohs(tcp_header.th_sport));
  net::AccessLayerPort dst_port(ntohs(tcp_header.th_dport));
  net::FiveTuple five_tuple(src_address, dst_address, net::kProtoTCP, src_port,
                            dst_port);

  uint16_t size = ntohs(ip_header.ip_len);
  uint16_t ip_id = ntohs(ip_header.ip_id);
  uint8_t ip_ttl = ip_header.ip_ttl;
  uint32_t seq = ntohl(tcp_header.th_seq);
  uint32_t ack = ntohl(tcp_header.th_ack);
  uint8_t tcp_flags = tcp_header.th_flags;

  UnidirectionalTCPFlowState& tcp_flow_state =
      tcp_flow_states.Emplace(five_tuple);

  tcp_flow_state.timestamps.Append(timestamp.count());
  tcp_flow_state.total_len_values.Append(size);
  tcp_flow_state.ip_id_values.Append(ip_id);
  tcp_flow_state.ttl_values.Append(ip_ttl);
  tcp_flow_state.sequence_numbers.Append(seq);
  tcp_flow_state.ack_numbers.Append(ack);
  tcp_flow_state.tcp_flag_values.Append(tcp_flags);
}

void FlowTracker::HandleUDP(Timestamp timestamp, const IPHeader& ip_header,
                            const UDPHeader& udp_header, uint16_t payload_len) {
  nc::Unused(timestamp);
  nc::Unused(ip_header);
  nc::Unused(udp_header);
  nc::Unused(payload_len);
}

void FlowTracker::HandleICMP(Timestamp timestamp, const IPHeader& ip_header,
                             const ICMPHeader& icmp_header,
                             uint16_t payload_len) {
  nc::Unused(timestamp);
  nc::Unused(ip_header);
  nc::Unused(icmp_header);
  nc::Unused(payload_len);
}

void FlowTracker::HandleUnknownIP(Timestamp timestamp,
                                  const IPHeader& ip_header,
                                  uint16_t payload_len) {
  nc::Unused(timestamp);
  nc::Unused(ip_header);
  nc::Unused(payload_len);
}

std::vector<BidirectionalTCPFlowState> FlowTracker::BidirectionalTCPFlows()
    const {
  std::unordered_map<FiveTuple, const UnidirectionalTCPFlowState*,
                     nc::net::FiveTupleHasher> values =
      tcp_flow_states.Values();

  std::vector<BidirectionalTCPFlowState> out;
  for (const auto& five_tuple_and_flow_state : values) {
    const FiveTuple& tuple = five_tuple_and_flow_state.first;
    const UnidirectionalTCPFlowState* flow_state =
        five_tuple_and_flow_state.second;

    FiveTuple reverse_tuple = tuple.Reverse();
    const UnidirectionalTCPFlowState* reverse_flow_state =
        nc::FindPtrOrNull(values, reverse_tuple);
    if (reverse_flow_state == nullptr) {
      continue;
    }

    if (reverse_tuple < tuple) {
      continue;
    }

    out.emplace_back(tuple, flow_state, reverse_flow_state);
  }

  return out;
}

void BidirectionalTCPFlowState::PopulateClientAndServer(
    const nc::net::FiveTuple& forward_tuple,
    const UnidirectionalTCPFlowState* forward_flow,
    const UnidirectionalTCPFlowState* reverse_flow) {
  // Will try to figure out which one is the server based on ports.
  if (forward_tuple.src_port() == nc::net::AccessLayerPort(80) ||
      forward_tuple.src_port() == nc::net::AccessLayerPort(443)) {
    client_to_server_tuple_ = forward_tuple.Reverse();
    client_to_server_flow_ = reverse_flow;
    server_to_client_flow_ = forward_flow;
    return;
  }

  client_to_server_tuple_ = forward_tuple;
  client_to_server_flow_ = forward_flow;
  server_to_client_flow_ = reverse_flow;
}

}  // namespace nc
