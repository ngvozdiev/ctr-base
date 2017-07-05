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

struct PcapDataTraceBin {
  PcapDataTraceBin() : bytes(0), packets(0), flows_enter(0), flows_exit(0) {}

  PcapDataTraceBin(const PBBin& bin_pb)
      : bytes(bin_pb.byte_count()),
        packets(bin_pb.packet_count()),
        flows_enter(bin_pb.enter_flow_count()),
        flows_exit(bin_pb.exit_flow_count()) {}

  PcapDataTraceBin(uint32_t bytes, uint32_t packets, uint32_t flows_enter,
                   uint32_t flows_exit)
      : bytes(bytes),
        packets(packets),
        flows_enter(flows_enter),
        flows_exit(flows_exit) {}

  // Combines this bin with a fraction of another.
  void Combine(const PcapDataTraceBin& other, double fraction);

  uint32_t bytes;
  uint32_t packets;
  uint32_t flows_enter;
  uint32_t flows_exit;
};

class PcapDataTrace;

// The index of a single slice from a trace.
struct TraceSliceTag {};
using TraceSliceIndex = nc::Index<TraceSliceTag, uint32_t>;
using TraceSliceSet = nc::PerfectHashSet<uint32_t, TraceSliceTag>;

// A bin sequence is a view into a slice of a trace (or multiple).
class BinSequence {
 public:
  struct TraceAndSlice {
    PcapDataTrace* trace;
    TraceSliceIndex slice;
    size_t start_bin;
    size_t end_bin;
    double fraction;
  };

  BinSequence(const std::vector<TraceAndSlice>& traces);

  BinSequence(const BinSequence& other) : BinSequence(other.traces_) {}

  // Combines this sequence's bins with another sequence's bins. Both
  // sequences should have the same bin size.
  void Combine(const BinSequence& other);

  // This sequence's bin size.
  const std::chrono::microseconds bin_size() const;

  size_t bin_count() const {
    size_t min = std::numeric_limits<size_t>::max();
    for (const auto& trace : traces_) {
      size_t bin_count = trace.end_bin - trace.start_bin;
      min = std::min(min, bin_count);
    }
    return min;
  }

  // Total bytes and packets in the whole sequence.
  std::pair<uint64_t, uint64_t> TotalBytesAndPackets();

  // Splits this sequence into a number of smaller sequences. Each will have the
  // same number of bins, but will only carry a part of this sequence's
  // sub-sequences. Fill die if there are more elements in 'fractions' than
  // there are sub-sequences. The sum of 'fractions' should be 1. Splitting is
  // deterministic and depends on the order of traces_.
  std::vector<BinSequence> SplitOrDie(
      const std::vector<double>& fractions) const;

  std::vector<BinSequence> PreciseSplitOrDie(
      const std::vector<double>& fractions) const;

  // Returns as many integers as there are bins. Given a rate, we can compute
  // how many bytes (B) will be transmitted for the period of a single bin. This
  // function returns B_i - B for each bin where B_i is the number of bytes in
  // the bin. If all values are positive this bin sequence will fit through a
  // link of rate 'rate' without causing congestion.
  std::vector<double> Residuals(nc::net::Bandwidth rate);

  // Generates an AggregateHistory from the entire range of this BinSequence.
  // Each of the history's bins will combine 'history_bin_size' bins from this
  // BinSequence.
  AggregateHistory GenerateHistory(std::chrono::milliseconds history_bin_size);

  // Returns another BinSequence with the same traces, but only contains the
  // first 'offset_from_start' bins.
  BinSequence CutFromStart(size_t offset_from_start) const;

  // Same as above, but takes in a duration. The duration should be an exact
  // multiple of bin_size.
  BinSequence CutFromStart(std::chrono::microseconds duration) const;

  // Offsets both the start and the end bin of this sequence.
  BinSequence Offset(size_t bin_count) const;

  // Optionally rebins the trace and returns the bins. If 'bin_size' is equal to
  // this trace's base bin size (returned by bin_size()) no rebinning will
  // happen. If it is a multiple of the base bin size will combine every N bins.
  std::vector<PcapDataTraceBin> AccumulateBins(
      std::chrono::microseconds bin_size);

  // Sum of all bytes (in bits) divided by bin_size * bin_count in seconds.
  nc::net::Bandwidth MeanRate();

 private:
  // Combines every 'bin_size_multiplier' bins into a PcapDataTraceBin and
  // returns the sequence.
  std::vector<PcapDataTraceBin> AccumulateBinsPrivate(
      size_t bin_size_multiplier);

  std::vector<TraceAndSlice> traces_;
};

// A collection of .pcap files that form a single trace.
class PcapDataTrace {
 public:
  static constexpr size_t kCacheSize = 10000000;
  static constexpr size_t kHintsSize = 1000;

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
  void Bins(TraceSliceIndex slice, size_t start_bin, size_t end_bin,
            std::function<void(const PBBin& binned_data)> callback);

  // Like Bins, but combines multiple slices into one. Will also pull the
  // results from the bin cache if possible.
  std::pair<std::vector<PcapDataTraceBin>::const_iterator,
            std::vector<PcapDataTraceBin>::const_iterator>
  BinsCombined(const TraceSliceSet& slices, size_t start_bin, size_t end_bin);

  TraceSliceSet AllSlices() const;

  // Constructs a BinSequence with slices from the trace.
  BinSequence ToSequence(const TraceSliceSet& slices);

  std::vector<std::string> trace_files() const;
  const TraceId& id() const { return id_; }
  const PBPcapDataTrace& ToProtobuf() const { return trace_pb_; }

  // All bins in this trace are by default this big.
  const std::chrono::microseconds base_bin_size() const;

  // Summarizes the trace in a human-readable form.
  std::string Summary();

  // Creates a new trace, at a given file at some offset.
  PcapDataTrace(const std::string& file, size_t offset);

  // Total size of the trace, its binned data and flows data in the file.
  size_t TotalSizeInFile() const;

 private:
  struct CachedBins {
    size_t from;
    size_t to;
    std::vector<PcapDataTraceBin> bins;
  };

  // Combines all bins for a set of slices and stores them in the cache.
  const CachedBins* AddToBinCache(const TraceSliceSet& slices, size_t from,
                                  size_t to);

  size_t TraceProtobufSize() const;

  // Locates this trace's entry in the store.
  std::string store_file_;
  size_t offset_in_file_;

  // A protobuf with the trace.
  PBPcapDataTrace trace_pb_;

  // The id.
  TraceId id_;

  // For each slice stores a bin number and an offset---bins with a bin number
  // larger or equal to this number are offset many bytes away from the start of
  // the slice's bins.
  std::map<TraceSliceIndex, std::map<size_t, size_t>> bin_start_hints_;

  // Cached, combined bins for series of slices.
  std::map<TraceSliceSet, CachedBins> bins_cache_;

  DISALLOW_COPY_AND_ASSIGN(PcapDataTrace);
};

// Stores and manages .pcap traces. Also caches expensive results.
class PcapTraceStore {
 public:
  PcapTraceStore(const std::string& file);

  // Gets a trace by id. The trace is owned by this object.
  PcapDataTrace* GetTraceOrNull(const TraceId& id);
  const PcapDataTrace* GetTraceOrNull(const TraceId& id) const;

  std::string Summary() const;

  std::set<TraceId> AllIds() const {
    std::set<TraceId> out;
    nc::InsertKeysFromMap(traces_, &out);
    return out;
  }

  std::vector<PcapDataTrace*> AllTraces() {
    std::vector<PcapDataTrace*> out;
    for (const auto& id_and_trace : traces_) {
      out.emplace_back(id_and_trace.second.get());
    }

    return out;
  }

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

// Will generate BinSequences given an initial list of BinSequences. If Next()
// is called more times than there are BinSequences in the initial list, will
// offset a previous BinSequence.
class BinSequenceGenerator {
 public:
  BinSequenceGenerator(std::vector<BinSequence>& initial_list,
                       size_t offset_step)
      : initial_list_(initial_list), offset_step_(offset_step), i_(0) {}

  BinSequence Next();

 private:
  std::vector<BinSequence> initial_list_;
  size_t offset_step_;
  size_t i_;
};

BinSequence BinsAtRate(nc::net::Bandwidth rate,
                       std::chrono::microseconds init_window,
                       BinSequenceGenerator* sequence_generator);

}  // namespace e2e

#endif
