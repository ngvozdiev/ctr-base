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
#include <string>
#include <tuple>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/packet.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_common.h"
#include "pcap_data.pb.h"
#include "flow_counter.h"

namespace ctr {
class AggregateHistory;
} /* namespace ctr */
namespace nc {
namespace htsim {
class PcapPacketGen;
} /* namespace htsim */
} /* namespace nc */

namespace ctr {

class TraceSubdivision {
 public:
  TraceSubdivision(size_t n, std::vector<bool> active_indices);
  TraceSubdivision(const PBTraceSubdivision& subdivision);

  // The original trace will be split N times.
  size_t n() const { return n_; }

  // Each flow of the original trace will end up in one of the N buckets, only
  // buckets with these indices will be considered.
  const std::vector<bool>& active_indices() const { return active_indices_; }

  friend bool operator<(const TraceSubdivision& a, const TraceSubdivision& b);

  std::string ToString() const;

  PBTraceSubdivision ToProtobuf() const;

 private:
  uint32_t n_;
  std::vector<bool> active_indices_;
};

// Uniquely identifies a trace.
class TraceId {
 public:
  TraceId(const std::string& location, size_t year, size_t month, size_t day,
          const std::vector<TraceSubdivision>& subdivisions = {});
  TraceId(const PBTraceId& trace_id);

  // Location where the trace was taken and time it was taken.
  const std::string& location() const { return location_; }

  std::tuple<uint16_t, uint16_t, uint16_t> timestamp() const {
    return std::make_tuple(year_, month_, day_);
  }

  const std::vector<TraceSubdivision>& subdivisions() const {
    return subdivisions_;
  }

  // Pretty-prints the trace id. This string uniquely identifies the trace, and
  // can be used as part of a filename (as long as the location string is safe).
  std::string ToString() const;

  std::string ToStringHumanReadable() const;

  // Pretty-prints the timestamp.
  std::string LocationAndTimestampToString() const;

  // Serializes this id to protobuf.
  PBTraceId ToProtobuf() const;

  // Creates a new TraceId for a trace that will split the original trace into N
  // chunks and will only see the chunks whose indices arein 'active_indices'.
  TraceId Subdivide(size_t n, const std::set<size_t>& active_indices) const;

  friend bool operator<(const TraceId& a, const TraceId& b);
  friend std::ostream& operator<<(std::ostream& output, const TraceId& op);

 private:
  // Location of the trace.
  std::string location_;

  // Time the trace was taken.
  uint16_t year_;
  uint16_t month_;
  uint16_t day_;

  // A list of subdivisions to be applied to the trace.
  std::vector<TraceSubdivision> subdivisions_;
};

class PcapTraceStore;

struct PcapDataTraceBin {
  uint64_t bytes = 0;
  uint64_t packets = 0;
  uint64_t flows = 0;
};

// A collection of .pcap files that form a single trace.
class PcapDataTrace {
 public:
  static constexpr size_t kCacheSize = 10000000;

  static std::vector<std::chrono::milliseconds> GetDefaultBinSizes();

  AggregateHistory InitialHistoryOrDie(
      std::chrono::milliseconds poll_period,
      std::chrono::milliseconds initial_window) const;

  // Records bins for all default bin sizes for the entire trace.
  void SinglePassCacheBins();

  // Records flow summaries for the entire trace.
  void SinglePassCacheFlows();

  // Bins the first 'duration' milliseconds from the trace and adds the bins to
  // the binned cache file.
  void BinCache(
      const std::vector<std::chrono::milliseconds>& bin_sizes,
      std::chrono::milliseconds duration = std::chrono::milliseconds::zero());

  // Returns num_bins bins from the start of the trace of a given size. If
  // the bins are not in the cache or will die.
  std::vector<PcapDataTraceBin> BinOrDie(std::chrono::milliseconds bin_size,
                                         size_t num_bins) const;

  std::vector<PcapDataTraceBin> BinOrEmpty(std::chrono::milliseconds bin_size,
                                           size_t num_bins) const;

  // Same as above, but returns all bins for a bin size in the cache.
  std::vector<PcapDataTraceBin> BinOrDie(
      std::chrono::milliseconds bin_size) const;

  // Returns a packet generator that can be used in simulation to get packets
  // from the trace.
  std::unique_ptr<nc::htsim::PcapPacketGen> GetPacketGenerator(
      nc::EventQueue* event_queue) const;

  // Returns the fraction of the volume of all TCP flows with SYN packets that
  // is due to flows that have an RTT of less than 'threshold'.
  double FractionOfVolumeLocal(std::chrono::milliseconds threshold);

  // Calls a function with all TCP flows in the trace with a SYN packet.
  void TCPSYNFlowsOrDie(
      std::function<void(const TCPSYNFlowSummary& syn_flow_summay)> callback)
      const;

  // Calls a function with all binned data in the cache.
  void AllBinsOrDie(
      std::function<void(const PBBinnedData& binned_data)> callback) const;

  // Creates a new trace or returns an existing trace from the store. The
  // new/old will split the original trace into N chunks and will only see the
  // chunks whose indices are in 'active_indices'. The returned trace will be
  // backed by the same pcap trace files.
  PcapDataTrace* Subdivide(size_t n, const std::set<size_t>& indices);
  const PcapDataTrace& SubdivideOrDie(size_t n,
                                      const std::set<size_t>& indices) const;

  std::vector<std::string> trace_files() const;
  const TraceId& id() const { return id_; }
  const PBPcapDataTrace& ToProtobuf() const { return trace_pb_; }

  // Updates the tag associated with this trace.
  void set_tag(const std::string& tag);
  const std::string& tag() const { return trace_pb_.tag(); }

  // If this function returns false then a call to TCPSYNFlows or
  // FractionOfVolumeLocal will have to cache all flows in the trace, which is
  // potentially very expensive.
  bool HasSynFlowsCached() const {
    return !trace_pb_.tcp_syn_flow_summary_file().empty();
  }

  // Prints a summary of this aggregate's bins.
  std::string BinsSummary() const;

 private:
  PcapDataTrace(const TraceId& id, const PBPcapDataTrace& trace_pb,
                PcapTraceStore* parent_store);

  TraceId id_;

  // A pointer to the parent store.
  PcapTraceStore* parent_store_;

  // A protobuf with the trace.
  PBPcapDataTrace trace_pb_;

  friend class PcapTraceStore;

  DISALLOW_COPY_AND_ASSIGN(PcapDataTrace);
};

// Stores and manages .pcap traces. Also caches expensive results.
class PcapTraceStore {
 public:
  PcapTraceStore(const std::string& file, bool ignore_missing = false);

  // Gets a trace by id. The trace is owned by this object.
  PcapDataTrace& GetTraceOrDie(const TraceId& id);
  const PcapDataTrace& GetTraceOrDie(const TraceId& id) const;

  // Adds a new trace to the store. If there is already a trace with the same id
  // and the same set of pcap_files will return that. If there is a trace with
  // the same id, but a different set of pcap_files will die.
  PcapDataTrace* AddTrace(const TraceId& id,
                          const std::vector<std::string>& pcap_files);

  // Returns the ids of all traces in the store whose tag matches a regex.
  std::vector<TraceId> AllIds(const std::string& tag_regex) const;

  // Writes out the contents of the store to the backing file. Will be called
  // automatically by all functions in the store or PcapDataTrace that modify
  // state.
  void Serialize();

  // Returns one or more traces whose max rates when combined will sum up to
  // roughly 'rate'. If rates closer to mean rates are needed the last argument
  // can be lowered closer to 0 (it it is 0 the mean rates will be used instead
  // of the max rates).
  std::vector<TraceId> IdsWithMaxRate(
      const std::vector<TraceId>& ids_to_consider, nc::net::Bandwidth rate,
      std::chrono::milliseconds poll_period,
      std::chrono::milliseconds init_window, double multiplier = 1.0);

 private:
  std::string file_;
  std::map<TraceId, std::unique_ptr<PcapDataTrace>> traces_;
};

// Bins a pcap trace.
class PcapTraceBinner : public nc::EventConsumer,
                        public nc::htsim::PacketObserver {
 public:
  PcapTraceBinner(std::chrono::milliseconds bin_size,
                  nc::EventQueue* event_queue)
      : nc::EventConsumer("binner", event_queue),
        bin_size_chrono_(bin_size),
        flow_counter_(event_queue) {
    bin_size_ = event_queue->ToTime(bin_size);
    CHECK(bin_size_ > nc::EventQueueTime::ZeroTime());
    EnqueueIn(bin_size_);
  }

  void HandleEvent() override {
    current_bin_.flows = flow_counter_.EstimateCount();
    bins_.emplace_back(current_bin_);

    current_bin_.bytes = 0;
    current_bin_.flows = 0;
    current_bin_.packets = 0;

    EnqueueIn(bin_size_);
  }

  void ObservePacket(const nc::htsim::Packet& pkt) override {
    current_bin_.bytes += pkt.size_bytes();
    current_bin_.packets += 1;
    flow_counter_.NewPacket(pkt.five_tuple());
  }

  const std::vector<PcapDataTraceBin>& bins() const { return bins_; }

  std::chrono::milliseconds bin_size() const { return bin_size_chrono_; }

 private:
  std::chrono::milliseconds bin_size_chrono_;
  PcapDataTraceBin current_bin_;
  std::vector<PcapDataTraceBin> bins_;
  nc::EventQueueTime bin_size_;
  FlowCounter flow_counter_;
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
