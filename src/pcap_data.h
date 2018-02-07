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

#include "ncode/event_queue.h"
#include "ncode/htsim/match.h"
#include "ncode/htsim/packet.h"
#include "ncode/htsim/pcap_consumer.h"
#include "ncode/net/net_common.h"
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

  uint32_t bytes;
  uint32_t packets;
  uint32_t flows_enter;
  uint32_t flows_exit;
};

// A variant of the struct above that only keeps track of bytes and flows
// entered.
struct __attribute__((packed)) TrimmedPcapDataTraceBin {
  TrimmedPcapDataTraceBin() : bytes(0), flows_enter(0) {}

  TrimmedPcapDataTraceBin(const PBBin& bin_pb)
      : bytes(bin_pb.byte_count()), flows_enter(bin_pb.enter_flow_count()) {
    CHECK(bin_pb.enter_flow_count() <= std::numeric_limits<uint32_t>::max());
  }

  void CombineWithFraction(const TrimmedPcapDataTraceBin& other,
                           double fraction);

  // Like above, but does not take a fraction.
  void Combine(const TrimmedPcapDataTraceBin& other);

  uint32_t bytes;
  uint32_t flows_enter;
};

class PcapDataTrace;
class PcapDataBinCache;

// The index of a single slice from a trace.
struct TraceSliceTag {};
using TraceSliceIndex = nc::Index<TraceSliceTag, uint32_t>;
using TraceSliceSet = nc::PerfectHashSet<uint32_t, TraceSliceTag>;

template <typename Value>
using TraceSliceMap = nc::PerfectHashMap<uint32_t, TraceSliceTag, Value>;

// A bin sequence is a view into a slice of a trace (or multiple).
class BinSequence {
 public:
  struct TraceAndSlice {
    TraceAndSlice(const PcapDataTrace* trace, TraceSliceIndex slice,
                  size_t start_bin, size_t end_bin, double precise_split)
        : trace(trace),
          slice(slice),
          start_bin(start_bin),
          end_bin(end_bin),
          precise_split(precise_split) {}

    const PcapDataTrace* trace;
    TraceSliceIndex slice;
    size_t start_bin;
    size_t end_bin;

    // Will assume all bins are from 'trace' * precise_split.
    double precise_split;
  };

  BinSequence(const std::vector<TraceAndSlice>& traces);

  BinSequence(std::vector<TraceAndSlice>::const_iterator from,
              std::vector<TraceAndSlice>::const_iterator to);

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
  std::pair<uint64_t, uint64_t> TotalBytesAndPackets() const;

  // Splits this sequence into a number of smaller sequences. Each will have the
  // same number of bins, but will only carry a part of this sequence's
  // sub-sequences.
  std::vector<std::unique_ptr<BinSequence>> PreciseSplitOrDie(
      const std::vector<double>& fractions) const;

  // Simluates how this sequence of bins would fit through a queue of given
  // rate. Returns the max queue size.
  std::chrono::milliseconds SimulateQueue(nc::net::Bandwidth rate,
                                          PcapDataBinCache* cache) const;

  // Generates an AggregateHistory from the entire range of this BinSequence.
  // Each of the history's bins will combine 'history_bin_size' bins from this
  // BinSequence.
  AggregateHistory GenerateHistory(std::chrono::milliseconds history_bin_size,
                                   uint64_t flow_count,
                                   PcapDataBinCache* cache) const;

  // Returns another BinSequence with the same traces, but only contains the
  // first 'offset_from_start' bins.
  std::unique_ptr<BinSequence> CutFromStart(size_t offset_from_start) const;

  // Same as above, but takes in a duration. The duration should be an exact
  // multiple of bin_size.
  std::unique_ptr<BinSequence> CutFromStart(
      std::chrono::microseconds duration) const;

  // Offsets both the start and the end bin of this sequence.
  std::unique_ptr<BinSequence> Offset(size_t bin_count) const;

  // Optionally rebins the trace and returns the bins. If 'bin_size' is equal to
  // this trace's base bin size (returned by bin_size()) no rebinning will
  // happen. If it is a multiple of the base bin size will combine every N bins.
  std::vector<TrimmedPcapDataTraceBin> AccumulateBins(
      std::chrono::microseconds bin_size, PcapDataBinCache* cache) const;

  // Sum of all bytes (in bits) divided by bin_size * bin_count in seconds.
  nc::net::Bandwidth MeanRate(PcapDataBinCache* cache) const;

  // The maximum instantaneous rate reached for any bin.
  nc::net::Bandwidth MaxRate(PcapDataBinCache* cache) const;

  // Serializes this BinSequence.
  PBBinSequence ToProtobuf() const;

  const std::vector<TraceAndSlice>& traces() const { return traces_; }

 private:
  // Combines every 'bin_size_multiplier' bins into a PcapDataTraceBin and
  // returns the sequence.
  std::vector<TrimmedPcapDataTraceBin> AccumulateBinsPrivate(
      size_t bin_size_multiplier, PcapDataBinCache* cache) const;

  std::vector<TraceAndSlice> traces_;

  DISALLOW_COPY_AND_ASSIGN(BinSequence);
};

// A collection of .pcap files that form a single trace.
class PcapDataTrace {
 public:
  static constexpr size_t kFlowCacheSize = 1000000;

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
  // (exclusive) of a single slice of the trace. Will read the bins from disk.
  void BinsFromDisk(
      TraceSliceIndex slice, size_t start_bin, size_t end_bin,
      std::function<void(const PBBin& binned_data)> callback) const;

  TraceSliceSet AllSlices() const;

  // Returns a vector of TraceAndSlice for a set of slices.
  std::vector<BinSequence::TraceAndSlice> TracesAndSlices(
      const TraceSliceSet& slices, std::chrono::milliseconds time_offset =
                                       std::chrono::milliseconds::zero()) const;

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

  // Size of trace_pb_.
  size_t trace_pb_size_;

  // The id.
  TraceId id_;

  DISALLOW_COPY_AND_ASSIGN(PcapDataTrace);
};

// Caches bins.
class PcapDataBinCache {
 public:
  PcapDataBinCache() {}

  // Returns the bins from a given slice.
  std::pair<std::vector<TrimmedPcapDataTraceBin>::const_iterator,
            std::vector<TrimmedPcapDataTraceBin>::const_iterator>
  Bins(const PcapDataTrace* trace, TraceSliceIndex slice, size_t start_bin,
       size_t end_bin);

  std::string CacheStats() const;

 private:
  // Number of bins in a cache line.
  static constexpr size_t kCacheLineBinCount = 100000;

  struct CachedTrace {
    size_t starting_bin = 0;
    std::vector<TrimmedPcapDataTraceBin> bins;
  };

  // The cached data.
  std::map<const PcapDataTrace*, TraceSliceMap<CachedTrace>> cached_bins_;

  DISALLOW_COPY_AND_ASSIGN(PcapDataBinCache);
};

// Stores and manages .pcap traces.
class PcapTraceStore {
 public:
  using FilterFunction = std::function<bool(const TraceId& id)>;

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

  std::vector<const PcapDataTrace*> AllTracesExcept(
      FilterFunction filter) const {
    std::vector<const PcapDataTrace*> out;
    for (const auto& id_and_trace : traces_) {
      const TraceId& id = id_and_trace.first;
      if (filter && filter(id)) {
        continue;
      }

      out.emplace_back(id_and_trace.second.get());
    }

    return out;
  }

  std::vector<const PcapDataTrace*> AllTraces() const {
    return AllTracesExcept(FilterFunction());
  }

  std::unique_ptr<BinSequence> BinSequenceFromProtobufOrDie(
      const PBBinSequence& sequence_pb) const;

  // Makes sure a each trace in a BinSequence extends until the end.
  std::unique_ptr<BinSequence> ExtendBinSequence(
      const BinSequence& bin_sequence) const;

 private:
  std::string file_;
  std::map<TraceId, std::unique_ptr<PcapDataTrace>> traces_;
};

// Stores information about what combination of traces fit a given rate.
class PcapTraceFitStore {
 public:
  static void AddToStore(const BinSequence& bin_sequence,
                         const std::string& output_file,
                         PcapDataBinCache* cache);

  PcapTraceFitStore(const std::string& file, const PcapTraceStore* store);

  // Returns the BinSequence that fits bw if none are found will return empty
  // unique pointer.
  std::unique_ptr<BinSequence> GetBinSequence(
      nc::net::Bandwidth bw,
      nc::net::Bandwidth threshold = nc::net::Bandwidth::FromBitsPerSecond(10));

 private:
  // For each rate one or more bin sequences that fit that rate.
  std::map<nc::net::Bandwidth, std::vector<std::unique_ptr<BinSequence>>>
      rate_to_bin_sequence_;

  const PcapTraceStore* store_;
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
  BinSequenceGenerator(const std::vector<const PcapDataTrace*>& all_traces,
                       const std::vector<std::chrono::milliseconds>& offsets);

  // Returns a BinSequence that will fit within the given 'target_rate' (will
  // have queues of 10ms or less) over 'init_window'.
  std::unique_ptr<BinSequence> Next(nc::net::Bandwidth target_rate,
                                    std::chrono::microseconds init_window,
                                    PcapDataBinCache* cache,
                                    std::mt19937* rnd) const;

 private:
  using TraceAndOffset =
      std::pair<const PcapDataTrace*, std::chrono::milliseconds>;

  std::vector<TraceAndOffset> all_traces_;
};

}  // namespace e2e

#endif
