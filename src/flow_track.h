#ifndef FLOW_TRACK_H
#define FLOW_TRACK_H

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
};

struct BidirectionalTCPFlowState {
  nc::pcap::Timestamp syn_seen_at;
  nc::pcap::Timestamp syn_ack_seen_at;
  nc::pcap::Timestamp fin_seen_at;
  nc::pcap::Timestamp reverse_fin_seen_at;

  const UnidirectionalTCPFlowState* forward_flow;
  const UnidirectionalTCPFlowState* reverse_flow;
};

class FlowTracker;

class UnidirectionalFlowStateCache
    : public nc::LRUCache<nc::net::FiveTuple, UnidirectionalFlowState,
                          nc::net::FiveTupleHasher> {
 public:
  UnidirectionalFlowStateCache(size_t cache_size, FlowTracker* parent)
      : nc::LRUCache<nc::net::FiveTuple, UnidirectionalFlowState,
                     nc::net::FiveTupleHasher>(cache_size),
        parent_(parent) {}

  void ItemEvicted(const nc::net::FiveTuple& key,
                   std::unique_ptr<UnidirectionalFlowState> value) override;

 private:
  FlowTracker* parent_;
};

class FlowTracker {
 private:
  UnidirectionalFlowStateCache flow_states;
};

}  // namespace nc

#endif
