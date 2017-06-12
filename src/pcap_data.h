#ifndef E2E_PCAP_DATA_H
#define E2E_PCAP_DATA_H

#include <stddef.h>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/htsim/packet.h"
#include "ncode_common/src/htsim/pcap_consumer.h"
#include "ncode_common/src/net/net_common.h"
#include "pcap_data.pb.h"

namespace ctr {
class AggregateHistory;
} /* namespace ctr */

namespace ctr {

// Uniquely identifies a trace.
class TraceId {
 public:
  TraceId(const std::string& location, size_t year, size_t month, size_t day);
  TraceId(const PBTraceId& trace_id);

  // Location where the trace was taken and time it was taken.
  const std::string& location() const { return location_; }

  std::tuple<uint16_t, uint16_t, uint16_t> timestamp() const {
    return std::make_tuple(year_, month_, day_);
  }

  // Pretty-prints the trace id.
  std::string ToString() const;

  // Serializes this id to protobuf.
  PBTraceId ToProtobuf() const;

  friend bool operator<(const TraceId& a, const TraceId& b);
  friend std::ostream& operator<<(std::ostream& output, const TraceId& op);

 private:
  // Location of the trace.
  std::string location_;

  // Time the trace was taken.
  uint16_t year_;
  uint16_t month_;
  uint16_t day_;
};

class PcapDataTraceBin {
 public:
  PcapDataTraceBin(const PBBin& bin_pb)
      : bytes_(bin_pb.byte_count()),
        packets_(bin_pb.packet_count()),
        flows_enter_(bin_pb.enter_flow_count()),
        flows_exit_(bin_pb.exit_flow_count()) {}

  PcapDataTraceBin(uint64_t bytes, uint64_t packets, uint64_t flows_enter,
                   uint64_t flows_exit)
      : bytes_(bytes),
        packets_(packets),
        flows_enter_(flows_enter),
        flows_exit_(flows_exit) {}

  uint64_t bytes() const { return bytes_; }

  uint64_t flows_enter() const { return flows_enter_; }

  uint64_t flows_exit() const { return flows_exit_; }

  uint64_t packets() const { return packets_; }

  // Combines this bin with another.
  void Combine(const PcapDataTraceBin& other);

 private:
  uint64_t bytes_ = 0;
  uint64_t packets_ = 0;
  uint64_t flows_enter_ = 0;
  uint64_t flows_exit_ = 0;
};

class PcapDataTrace;

// A bin sequence is a view into a slice of a trace (or multiple).
class BinSequence {
 public:
  using TraceAndSlice = std::pair<const PcapDataTrace*, size_t>;
  BinSequence(const std::vector<TraceAndSlice>& traces, size_t start_bin,
              size_t end_bin);

  // Combines this sequence's bins with another sequence's bins. Both
  // sequences should have the same number of bins and bin size.
  void Combine(const BinSequence& other);

  // This sequence's bin size.
  const std::chrono::microseconds bin_size() const;

  size_t bin_count() const;

  // Total bytes and packets in the whole sequence.
  std::pair<uint64_t, uint64_t> TotalBytesAndPackets() const;

  // Splits this sequence into a number of smaller sequences. Each will have the
  // same number of bins, but will only carry a part of this sequence's
  // sub-sequences. Fill die if there are more elements in 'fractions' than
  // there are sub-sequences. The sum of 'fractions' should be 1. Splitting is
  // deterministic and depends on the order of traces_.
  std::vector<BinSequence> SplitOrDie(
      const std::vector<double>& fractions) const;

  // Returns as many integers as there are bins. Given a rate, we can compute
  // how many bytes (B) will be transmitted for the period of a single bin. This
  // function returns B_i - B for each bin where B_i is the number of bytes in
  // the bin. If all values are positive this bin sequence will fit through a
  // link of rate 'rate' without causing congestion.
  std::vector<double> Residuals(nc::net::Bandwidth rate) const;

  // Generates an AggregateHistory from the entire range of this BinSequence.
  // Each of the history's bins will combine 'history_bin_size' bins from this
  // BinSequence.
  AggregateHistory GenerateHistory(
      std::chrono::milliseconds history_bin_size) const;

  // Returns another BinSequence with the same traces, but a smaller range.
  BinSequence LimitRange(size_t start_bin, size_t end_bin) const;

 private:
  std::vector<PcapDataTraceBin> AccumulateBins() const;

  // Defines a range.
  size_t start_bin_;
  size_t end_bin_;

  std::vector<TraceAndSlice> traces_;
};

// A collection of .pcap files that form a single trace.
class PcapDataTrace {
 public:
  static constexpr size_t kCacheSize = 10000000;

  // Constructs a PcapDataTrace from a partially populated 'trace_pb'. Will
  // append persistent information to 'output_file'.
  static void Init(const PBPcapDataTrace& trace_pb,
                   const std::string& output_file);

  // Returns a packet generator that can be used in simulation to get packets
  // from the trace.
  std::unique_ptr<nc::htsim::PcapPacketGen> GetPacketGenerator(
      nc::EventQueue* event_queue) const;

  // Calls a function with all TCP flows in the trace with a SYN packet.
  void TCPSYNFlows(std::function<void(const TCPSYNFlowSummary& syn_flow_summay)>
                       callback) const;

  // Calls a function with bins between stat_bin (inclusive) and end_bin
  // (exclusive) of a single slice of the trace.
  void Bins(size_t slice, size_t start_bin, size_t end_bin,
            std::function<void(const PBBin& binned_data)> callback) const;

  std::set<size_t> AllSlices() const;

  // Constructs a BinSequence with slices from the trace.
  BinSequence ToSequence(const std::set<size_t>& slices) const;

  std::vector<std::string> trace_files() const;
  const TraceId& id() const { return id_; }
  const PBPcapDataTrace& ToProtobuf() const { return trace_pb_; }

  // All bins in this trace are by default this big.
  const std::chrono::microseconds base_bin_size() const;

  // Summarizes the trace in a human-readable form.
  std::string Summary() const;

  // Creates a new trace, at a given file at some offset.
  PcapDataTrace(const std::string& file, size_t offset);

  // Total size of the trace, its binned data and flows data in the file.
  size_t TotalSizeInFile() const;

 private:
  size_t TraceProtobufSize() const;

  // Locates this trace's entry in the store.
  std::string store_file_;
  size_t offset_in_file_;

  // A protobuf with the trace.
  PBPcapDataTrace trace_pb_;

  // The id.
  TraceId id_;

  DISALLOW_COPY_AND_ASSIGN(PcapDataTrace);
};

// Stores and manages .pcap traces. Also caches expensive results.
class PcapTraceStore {
 public:
  PcapTraceStore(const std::string& file);

  // Gets a trace by id. The trace is owned by this object.
  PcapDataTrace& GetTraceOrDie(const TraceId& id);
  const PcapDataTrace& GetTraceOrDie(const TraceId& id) const;

  std::string Summary() const;

 private:
  std::string file_;
  std::map<TraceId, std::unique_ptr<PcapDataTrace>> traces_;
};

// Bins a pcap trace.
class PcapTraceBinner : public nc::EventConsumer,
                        public nc::htsim::PacketObserver {
 public:
  PcapTraceBinner(std::chrono::microseconds bin_size,
                  nc::EventQueue* event_queue)
      : nc::EventConsumer("binner", event_queue),
        bytes_(0),
        pkts_(0),
        enter_(0),
        exit_(0),
        bin_size_chrono_(bin_size) {
    bin_size_ = event_queue->ToTime(bin_size);
    CHECK(bin_size_ > nc::EventQueueTime::ZeroTime());
    EnqueueIn(bin_size_);
  }

  void HandleEvent() override {
    bins_.emplace_back(bytes_, pkts_, enter_, exit_);
    bytes_ = 0;
    pkts_ = 0;
    enter_ = 0;
    exit_ = 0;

    EnqueueIn(bin_size_);
  }

  void ObservePacket(const nc::htsim::Packet& pkt) override {
    bytes_ += pkt.size_bytes();
    pkts_ += 1;

    if (pkt.is_tcp()) {
      const nc::htsim::TCPPacket* tcp_packet =
          static_cast<const nc::htsim::TCPPacket*>(&pkt);
      if (tcp_packet->syn()) {
        enter_ += 1;
      } else if (tcp_packet->fin()) {
        exit_ += 1;
      }
    }
  }

  const std::vector<PcapDataTraceBin>& bins() const { return bins_; }

  std::chrono::microseconds bin_size() const { return bin_size_chrono_; }

 private:
  uint64_t bytes_;
  uint64_t pkts_;
  uint64_t enter_;
  uint64_t exit_;

  std::chrono::microseconds bin_size_chrono_;
  std::vector<PcapDataTraceBin> bins_;
  nc::EventQueueTime bin_size_;
};

// Consistently hashes incoming packets and evenly spreads them betwenn N
// observers.
class BalancingObserver : public nc::htsim::PacketObserver {
 public:
  BalancingObserver(const std::vector<nc::htsim::PacketObserver*>& observers);

  void ObservePacket(const nc::htsim::Packet& pkt) override;

 private:
  std::unique_ptr<nc::htsim::MatchRule> dispatch_rule_;
  std::vector<nc::htsim::PacketObserver*> observers_;
};

// A series of observers. Each packet is passed through all of them.
class ObserverPack : public nc::htsim::PacketHandler {
 public:
  ObserverPack(const std::vector<nc::htsim::PacketObserver*>& observers)
      : observers_(observers) {}

  void HandlePacket(nc::htsim::PacketPtr pkt) {
    for (nc::htsim::PacketObserver* observer : observers_) {
      observer->ObservePacket(*pkt);
    }
  }

 private:
  std::vector<nc::htsim::PacketObserver*> observers_;
};

}  // namespace e2e

#endif
