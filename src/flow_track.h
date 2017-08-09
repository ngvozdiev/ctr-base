#ifndef FLOW_TRACK_H
#define FLOW_TRACK_H

#include <chrono>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include "ncode_common/src/lru_cache.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/pcap.h"
#include "ncode_common/src/packer.h"

namespace nc {

// A TCP or UDP unidirectional flow.
struct UnidirectionalFlowState {
  nc::PackedUintSeq timestamps;
  nc::RLEField<uint16_t> total_len_values;
  nc::RLEField<uint8_t> ttl_values;
  nc::RLEField<uint16_t> ip_id_values;
};

struct UnidirectionalTCPFlowState : public UnidirectionalFlowState {
  nc::RLEField<uint32_t> sequence_numbers;
  nc::RLEField<uint32_t> ack_numbers;
  nc::RLEField<uint8_t> tcp_flag_values;
};

struct DataCluster {
  DataCluster(nc::pcap::Timestamp first_byte_at)
      : first_byte_at(first_byte_at),
        last_byte_at(first_byte_at),
        packets(0),
        bytes(0) {}

  std::string ToString() const;

  nc::pcap::Timestamp first_byte_at;
  nc::pcap::Timestamp last_byte_at;
  uint64_t packets;
  uint64_t bytes;
};

class BidirectionalTCPFlowState {
 public:
  BidirectionalTCPFlowState(const nc::net::FiveTuple& forward_tuple,
                            const UnidirectionalTCPFlowState* forward_flow,
                            const UnidirectionalTCPFlowState* reverse_flow) {
    PopulateClientAndServer(forward_tuple, forward_flow, reverse_flow);
  }

  // Splits the TCP flow into bursts of data that are at least max_delta apart.
  std::vector<DataCluster> BreakDown(nc::pcap::Timestamp max_delta) const;

  // Estimates the RTT of the connection from the initial SYN/SYN+ACK handshake.
  // If the handshake is not captured will return Timestamp::max().
  nc::pcap::Timestamp EstimateRTT() const;

  std::string ToString() const;

 private:
  void PopulateClientAndServer(const nc::net::FiveTuple& forward_tuple,
                               const UnidirectionalTCPFlowState* forward_flow,
                               const UnidirectionalTCPFlowState* reverse_flow);

  nc::net::FiveTuple client_to_server_tuple_;
  const UnidirectionalTCPFlowState* client_to_server_flow_;
  const UnidirectionalTCPFlowState* server_to_client_flow_;
};

class FlowTracker : public nc::pcap::PacketHandler {
 public:
  FlowTracker(size_t max_flow_count) : tcp_flow_states(max_flow_count) {}

  void HandleTCP(nc::pcap::Timestamp timestamp,
                 const nc::pcap::IPHeader& ip_header,
                 const nc::pcap::TCPHeader& tcp_header,
                 uint16_t payload_len) override;

  void HandleUDP(nc::pcap::Timestamp timestamp,
                 const nc::pcap::IPHeader& ip_header,
                 const nc::pcap::UDPHeader& udp_header,
                 uint16_t payload_len) override;

  void HandleICMP(nc::pcap::Timestamp timestamp,
                  const nc::pcap::IPHeader& ip_header,
                  const nc::pcap::ICMPHeader& icmp_header,
                  uint16_t payload_len) override;

  void HandleUnknownIP(nc::pcap::Timestamp timestamp,
                       const nc::pcap::IPHeader& ip_header,
                       uint16_t payload_len) override;

  std::vector<BidirectionalTCPFlowState> BidirectionalTCPFlows() const;

 private:
  nc::LRUCache<nc::net::FiveTuple, UnidirectionalTCPFlowState,
               nc::net::FiveTupleHasher> tcp_flow_states;
};

}  // namespace nc

#endif
