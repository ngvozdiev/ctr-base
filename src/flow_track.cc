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
using namespace std::chrono;

using TimeAndSize = std::pair<uint64_t, uint16_t>;

std::string DataCluster::ToString() const {
  nc::pcap::Timestamp duration = last_byte_at - first_byte_at;
  uint64_t duration_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
  return nc::StrCat("duration: ", duration_ms, "ms, pkts: ", packets,
                    " bytes: ", bytes);
}

static std::vector<TimeAndSize> GetTimesAndSizes(
    const UnidirectionalTCPFlowState& flow_state) {
  std::vector<uint64_t> timestamps_forward = flow_state.timestamps.Restore();
  std::vector<uint16_t> sizes_forward = flow_state.total_len_values.Restore();

  std::vector<TimeAndSize> out;
  out.reserve(timestamps_forward.size());

  CHECK(timestamps_forward.size() == sizes_forward.size());
  for (size_t i = 0; i < timestamps_forward.size(); ++i) {
    out.emplace_back(timestamps_forward[i], sizes_forward[i]);
  }

  return out;
}

static std::vector<TimeAndSize> Combine(
    const UnidirectionalTCPFlowState& forward,
    const UnidirectionalTCPFlowState& reverse) {
  std::vector<TimeAndSize> data_one = GetTimesAndSizes(forward);
  std::vector<TimeAndSize> data_two = GetTimesAndSizes(reverse);

  std::vector<TimeAndSize> all_data;
  std::move(data_one.begin(), data_one.end(), std::back_inserter(all_data));
  std::move(data_two.begin(), data_two.end(), std::back_inserter(all_data));

  std::sort(all_data.begin(), all_data.end(),
            [](const TimeAndSize& lhs, const TimeAndSize& rhs) {
              return lhs.first < rhs.first;
            });
  return all_data;
}

static std::vector<DataCluster> BreakDownStatic(
    Timestamp max_delta, const UnidirectionalTCPFlowState& forward,
    const UnidirectionalTCPFlowState& reverse) {
  std::vector<TimeAndSize> data_combined = Combine(forward, reverse);

  std::vector<DataCluster> clusters;
  for (size_t i = 0; i < data_combined.size(); ++i) {
    TimeAndSize time_and_size = data_combined[i];
    uint64_t time = time_and_size.first;
    uint16_t size = time_and_size.second;

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
    current_cluster.bytes += size;
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
  Timestamp rtt_estimate = EstimateRTT();
  std::vector<std::string> cluster_strings;
  for (const DataCluster& cluster : clusters) {
    cluster_strings.emplace_back(cluster.ToString());
  }

  return StrCat(out, Join(cluster_strings, ", "), "\nRTT estimate: ",
                duration_cast<milliseconds>(rtt_estimate).count(), "ms");
}

Timestamp BidirectionalTCPFlowState::EstimateRTT() const {
  std::vector<uint8_t> flags =
      client_to_server_flow_->tcp_flag_values.Restore();
  CHECK(!flags.empty());
  uint8_t flags_of_first_packet = flags.front();
  if (!(flags_of_first_packet & TCPHeader::kSynFlag)) {
    return nc::pcap::Timestamp::max();
  }

  std::vector<uint8_t> reverse_flags =
      server_to_client_flow_->tcp_flag_values.Restore();
  CHECK(!reverse_flags.empty());
  uint8_t flags_of_first_reverse_packet = reverse_flags.front();
  if (!((flags_of_first_reverse_packet & TCPHeader::kSynFlag) &&
        (flags_of_first_reverse_packet & TCPHeader::kAckFlag))) {
    return nc::pcap::Timestamp::max();
  }

  uint64_t time_of_first_packet =
      client_to_server_flow_->timestamps.Restore().front();
  uint64_t time_of_first_reverse_packet =
      server_to_client_flow_->timestamps.Restore().front();
  CHECK(time_of_first_reverse_packet > time_of_first_packet);
  return Timestamp(time_of_first_reverse_packet - time_of_first_packet);
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
