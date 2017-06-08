#include "pcap_data.h"

#include <fcntl.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <sys/stat.h>
#include <algorithm>
#include <cerrno>
#include <cstring>
#include <iterator>
#include <regex>
#include <utility>

#include "ncode_common/src/file.h"
#include "ncode_common/src/htsim/bulk_gen.h"
#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/htsim/pcap_consumer.h"
#include "ncode_common/src/lru_cache.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/md5.h"
#include "ncode_common/src/net/pcap.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/substitute.h"
#include "common.h"

namespace ctr {

struct TCPFlowRecord {
  TCPFlowRecord(size_t bytes, uint8_t ttl, bool forward)
      : forward(forward),
        rtt(0),
        bytes_seen(bytes),
        pkts_seen(1),
        ttl(ttl),
        rtt_fin(0) {}

  // True if the first packet from this flow that we saw was a SYN (as opposed
  // to a SYN/ACK).
  bool forward;
  std::chrono::nanoseconds rtt;
  size_t bytes_seen;
  size_t pkts_seen;
  uint8_t ttl;

  // An RTT estimate based on the FIN. Not always present and set to 0 if not
  // present.
  std::chrono::nanoseconds rtt_fin;
};

struct TCPFlowState {
  TCPFlowState(nc::EventQueueTime syn_seen_at, size_t bytes, uint8_t ttl,
               bool forward)
      : syn_seen_at(syn_seen_at),
        fin_seen_at(nc::EventQueueTime::ZeroTime()),
        fin_seq(0),
        flow_record(bytes, ttl, forward) {}

  nc::EventQueueTime syn_seen_at;
  nc::EventQueueTime fin_seen_at;

  // If we see a FIN for a connection we may after it see more ACKs for other
  // data, need to disambiguate in order to get the right RTT.
  uint32_t fin_seq;

  TCPFlowRecord flow_record;
};

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

static std::unique_ptr<google::protobuf::io::FileOutputStream> GetOutputStream(
    const std::string& filename, bool append) {
  int fd;

  if (append) {
    fd = open(filename.c_str(), O_WRONLY | O_APPEND | O_CREAT,
              S_IREAD | S_IWRITE | S_IRGRP | S_IROTH | S_ISUID);
  } else {
    fd = open(filename.c_str(), O_WRONLY | O_TRUNC | O_CREAT,
              S_IREAD | S_IWRITE | S_IRGRP | S_IROTH | S_ISUID);
  }

  CHECK(fd != -1) << "Bad output file " << filename;
  return nc::make_unique<google::protobuf::io::FileOutputStream>(fd);
}

static std::unique_ptr<google::protobuf::io::FileInputStream> GetInputStream(
    const std::string& filename) {
  int fd = open(filename.c_str(), O_RDONLY);
  CHECK(fd > 0) << "Bad input file " << filename << ": " << strerror(errno);
  return nc::make_unique<google::protobuf::io::FileInputStream>(fd);
}

class TCPFlowStateCache : public nc::LRUCache<nc::net::FiveTuple, TCPFlowState,
                                              nc::net::FiveTupleHasher> {
 public:
  ~TCPFlowStateCache() {
    // close streams
    file_output_->Close();
    file_output_.reset();
  }

  TCPFlowStateCache(const std::string& output_file, size_t cache_size)
      : nc::LRUCache<nc::net::FiveTuple, TCPFlowState,
                     nc::net::FiveTupleHasher>(cache_size),
        count_(0) {
    file_output_ = GetOutputStream(output_file, false);
  }

  void ItemEvicted(const nc::net::FiveTuple& key,
                   std::unique_ptr<TCPFlowState> value) override {
    const TCPFlowRecord& flow_record = value->flow_record;
    if (flow_record.pkts_seen == 1) {
      return;
    }

    if (flow_record.rtt < std::chrono::microseconds(1)) {
      return;
    }

    TCPSYNFlowSummary summary;
    summary.set_bytes_seen(flow_record.bytes_seen);
    summary.set_pkts_seen(flow_record.pkts_seen);
    summary.set_destination(key.ip_dst().Raw());
    summary.set_source(key.ip_src().Raw());
    summary.set_src_port(key.src_port().Raw());
    summary.set_dst_port(key.dst_port().Raw());
    summary.set_first_packet_ttl(flow_record.ttl);
    summary.set_pure_syn_seen(flow_record.forward);
    summary.set_syn_rtt_estimate_nanos(flow_record.rtt.count());
    summary.set_fin_rtt_estimate_nanos(flow_record.rtt_fin.count());
    CHECK(WriteDelimitedTo(summary, file_output_.get()));

    ++count_;
  }

  size_t count() const { return count_; }

 private:
  size_t count_;

  // File output stream. Owned by this object.
  std::unique_ptr<google::protobuf::io::FileOutputStream> file_output_;
};

class TCPFlowStateRecorder : public nc::htsim::PacketObserver {
 public:
  TCPFlowStateRecorder(const std::string& output_file, size_t cache_size,
                       nc::EventQueue* event_queue)
      : total_bytes_(0),
        total_pkts_(0),
        cache_(output_file, cache_size),
        event_queue_(event_queue) {}

  void ObservePacket(const nc::htsim::Packet& pkt) override {
    total_bytes_ += pkt.size_bytes();
    ++total_pkts_;

    if (pkt.five_tuple().ip_proto() == nc::net::kProtoTCP) {
      const nc::htsim::TCPPacket* tcp_packet =
          static_cast<const nc::htsim::TCPPacket*>(&pkt);

      bool syn = (tcp_packet->flags() & nc::pcap::TCPHeader::kSynFlag);
      bool ack = (tcp_packet->flags() & nc::pcap::TCPHeader::kAckFlag);
      bool fin = (tcp_packet->flags() & nc::pcap::TCPHeader::kFinFlag);
      AddPacket(pkt, tcp_packet->sequence().Raw(), syn, syn && !ack, fin);
    }
  }

  size_t FinishUp() {
    cache_.EvictAll();
    return cache_.count();
  }

  size_t total_bytes() const { return total_bytes_; }

  size_t total_pkts() const { return total_pkts_; }

 private:
  void AddPacket(const nc::htsim::Packet& packet, uint32_t seq, bool first_pkt,
                 bool forward, bool fin) {
    uint8_t ttl = packet.ttl();
    uint16_t bytes = packet.size_bytes();
    const nc::net::FiveTuple& five_tuple = packet.five_tuple();

    nc::EventQueueTime now = event_queue_->CurrentTime();
    if (first_pkt) {
      // If there is already state for the same 5-tuple then we will replace it.
      cache_.InsertNew(five_tuple, now, bytes, ttl, forward);
      return;
    }

    TCPFlowState* flow_state_ptr = cache_.FindOrNull(five_tuple);
    if (flow_state_ptr == nullptr) {
      // No SYN recorded for this flow.
      return;
    }

    TCPFlowRecord& flow_record = flow_state_ptr->flow_record;
    if (flow_record.pkts_seen == 1) {
      // This is the second packet from the flow -- the one after the syn. Will
      // update the RTT estimate.
      nc::EventQueueTime time_syn_seen = flow_state_ptr->syn_seen_at;
      CHECK(now >= time_syn_seen);
      nc::EventQueueTime rtt_estimate = now - time_syn_seen;
      auto rtt_est = event_queue_->TimeToNanos(rtt_estimate);
      if (rtt_est < std::chrono::microseconds(1)) {
        LOG(ERROR) << "Sub-microsecond RTT estimate, pkt: " << packet.ToString()
                   << " RTT est " << rtt_est.count() << " nanos";
      }
      flow_record.rtt = rtt_est;
    }

    nc::EventQueueTime time_fin_seen = flow_state_ptr->fin_seen_at;
    if (fin) {
      if (time_fin_seen != nc::EventQueueTime::ZeroTime()) {
        if (flow_state_ptr->fin_seq == seq) {
          // It is a retransmission. The RTT estimate will be wrong. Better to
          // reset it to 0 so that we ignore it.
          flow_state_ptr->fin_seq = 0;
          flow_record.rtt_fin = std::chrono::nanoseconds(0);
        }
      } else {
        flow_state_ptr->fin_seen_at = now;
        flow_state_ptr->fin_seq = seq;
      }
    } else {
      if (time_fin_seen != nc::EventQueueTime::ZeroTime() &&
          seq == flow_state_ptr->fin_seq + 1) {
        CHECK(now >= time_fin_seen);
        nc::EventQueueTime rtt_estimate = now - time_fin_seen;
        flow_record.rtt_fin = event_queue_->TimeToNanos(rtt_estimate);
      }
    }

    ++flow_record.pkts_seen;
    flow_record.bytes_seen += bytes;
  }

  // Total number of bytes/pkts.
  size_t total_bytes_;
  size_t total_pkts_;

  // Contains information about flows.
  TCPFlowStateCache cache_;

  // The main event queue.
  nc::EventQueue* event_queue_;
};

static std::string ExtractDirName(const std::string& file_location) {
  std::vector<std::string> pieces = nc::Split(file_location, "/", true);
  CHECK(pieces.size() > 0);

  std::string out;
  nc::Join(pieces.begin(), std::prev(pieces.end()), "/", &out);
  return nc::StrCat("/", out);
}

template <typename T>
static bool ReadDelimitedFrom(
    T* message, google::protobuf::io::FileInputStream* file_input) {
  // We create a new coded stream for each message.
  google::protobuf::io::CodedInputStream input(file_input);

  // Read the size.
  uint32_t size;
  if (!input.ReadVarint32(&size)) {
    return false;
  }

  // Tell the stream not to read beyond that size.
  google::protobuf::io::CodedInputStream::Limit limit = input.PushLimit(size);

  // Parse the message.
  if (!message->MergeFromCodedStream(&input)) {
    return false;
  }

  if (!input.ConsumedEntireMessage()) {
    return false;
  }

  // Release the limit.
  input.PopLimit(limit);

  return true;
}

static std::unique_ptr<nc::pcap::OfflineSourceProvider> GetProvider(
    const PcapDataTrace& data_trace) {
  std::vector<std::string> files;
  for (const std::string& file : data_trace.trace_files()) {
    files.emplace_back(file);
  }

  std::unique_ptr<nc::pcap::OfflineSourceProvider> offline_provider =
      nc::make_unique<nc::pcap::DefaultOfflineSourceProvider>(files);
  return offline_provider;
}

class DownscalingPcapPacketGen : public nc::htsim::PcapPacketGen {
 public:
  DownscalingPcapPacketGen(const PcapDataTrace& pcap_data_trace,
                           nc::EventQueue* event_queue)
      : nc::htsim::PcapPacketGen(GetProvider(pcap_data_trace), event_queue) {
    using namespace nc::htsim;

    const TraceId& id = pcap_data_trace.id();

    for (const TraceSubdivision& subdivision : id.subdivisions()) {
      nc::htsim::MatchRuleKey key(kWildPacketTag, kWildDevicePortNumber,
                                  {nc::net::FiveTuple::kDefaultTuple});
      auto rule = nc::make_unique<MatchRule>(key);

      // Will add N actions to the rule, all with the same weight, this will
      // split incoming traffic N ways.
      for (size_t i = 0; i < subdivision.n(); ++i) {
        auto action = nc::make_unique<MatchRuleAction>(
            nc::net::DevicePortNumber(i), kWildPacketTag, 1);
        rule->AddAction(std::move(action));
      }

      DownscaleState downscale_state;
      downscale_state.downscale_rule = std::move(rule);
      downscale_state.ports = subdivision.active_indices();
      downscale_states_.emplace_back(std::move(downscale_state));
    }

    set_max_inter_packet_gap(std::chrono::duration_cast<nc::pcap::Timestamp>(
        std::chrono::milliseconds(10)));
  }

  bool Ignore(const nc::net::FiveTuple& five_tuple) override {
    for (DownscaleState& downscale_state : downscale_states_) {
      nc::htsim::MatchRuleAction* action =
          downscale_state.downscale_rule->ChooseOrNull(five_tuple);
      CHECK(action != nullptr);

      uint32_t index = action->output_port().Raw();
      CHECK(index < downscale_state.ports.size());
      if (!downscale_state.ports[index]) {
        return true;
      }
    }

    return false;
  }

 private:
  struct DownscaleState {
    std::unique_ptr<nc::htsim::MatchRule> downscale_rule;
    std::vector<bool> ports;
  };

  std::vector<DownscaleState> downscale_states_;
};

std::vector<std::chrono::milliseconds> PcapDataTrace::GetDefaultBinSizes() {
  using namespace std::chrono;
  return {milliseconds(10), milliseconds(100), milliseconds(1000),
          milliseconds(10000), milliseconds(60000)};
}

PcapDataTrace::PcapDataTrace(const TraceId& id, const PBPcapDataTrace& trace_pb,
                             PcapTraceStore* parent_store)
    : id_(id), parent_store_(parent_store), trace_pb_(trace_pb) {
  *trace_pb_.mutable_id() = id.ToProtobuf();
}

std::unique_ptr<nc::htsim::PcapPacketGen> PcapDataTrace::GetPacketGenerator(
    nc::EventQueue* event_queue) const {
  std::unique_ptr<nc::htsim::PcapPacketGen> pcap_packet_gen =
      nc::make_unique<DownscalingPcapPacketGen>(*this, event_queue);
  return pcap_packet_gen;
}

void PcapDataTrace::TCPSYNFlowsOrDie(
    std::function<void(const TCPSYNFlowSummary& syn_flow_summay)> callback)
    const {
  // Will check if there is a file with the TCP flows.
  std::string tcp_flows_file = trace_pb_.tcp_syn_flow_summary_file();
  CHECK(!tcp_flows_file.empty());
  if (!callback) {
    return;
  }

  TCPSYNFlowSummary flow_summary;
  auto file_input = GetInputStream(tcp_flows_file);
  while (ReadDelimitedFrom(&flow_summary, file_input.get())) {
    callback(flow_summary);
    flow_summary.Clear();
  }
}

void PcapDataTrace::AllBinsOrDie(
    std::function<void(const PBBinnedData& binned_data)> callback) const {
  std::string binned_data_file = trace_pb_.binned_data_file();
  CHECK(!binned_data_file.empty());
  if (!callback) {
    return;
  }

  PBBinnedData binned_data_pb;
  auto file_input = GetInputStream(binned_data_file);
  while (ReadDelimitedFrom(&binned_data_pb, file_input.get())) {
    callback(binned_data_pb);
    binned_data_pb.Clear();
  }
  file_input->Close();
}

std::string PcapDataTrace::BinsSummary() const {
  std::string out = nc::StrCat(id_.ToStringHumanReadable(), ":\n");
  AllBinsOrDie([&out](const PBBinnedData& binned_data) {
    uint64_t bin_size_ms = binned_data.bin_size_ms();
    uint64_t bin_count = binned_data.bytes_binned_size();
    nc::StrAppend(&out, "\tbin ", bin_size_ms, "ms (", bin_count, " bins)\n");
  });

  return out;
}

PcapDataTrace* PcapDataTrace::Subdivide(
    size_t n, const std::set<size_t>& active_indices) {
  TraceId new_id = id_.Subdivide(n, active_indices);
  return parent_store_->AddTrace(new_id, trace_files());
}

const PcapDataTrace& PcapDataTrace::SubdivideOrDie(
    size_t n, const std::set<size_t>& active_indices) const {
  TraceId new_id = id_.Subdivide(n, active_indices);
  return parent_store_->GetTraceOrDie(new_id);
}

AggregateHistory PcapDataTrace::InitialHistoryOrDie(
    std::chrono::milliseconds poll_period,
    std::chrono::milliseconds initial_window) const {
  uint64_t poll_period_ms = poll_period.count();
  uint64_t initial_window_ms = initial_window.count();
  size_t bins_count = initial_window_ms / poll_period_ms;

  // Will get the first N bins.
  std::vector<PcapDataTraceBin> all_bins = BinOrDie(poll_period, bins_count);
  CHECK(all_bins.size() == bins_count);
  std::vector<uint64_t> byte_bins;
  for (size_t i = 0; i < bins_count; ++i) {
    byte_bins.emplace_back(all_bins[i].bytes);
  }

  // Will get the number of flows in FLAGS_pcap_trace_flow_estimate_window_ms
  std::vector<PcapDataTraceBin> all_periods = BinOrDie(initial_window, 1);
  CHECK(!all_periods.empty());
  uint64_t flow_count = all_periods.front().flows;

  return AggregateHistory(byte_bins, poll_period, flow_count);
}

std::vector<PcapDataTraceBin> PcapDataTrace::BinOrDie(
    std::chrono::milliseconds bin_size, size_t num_bins) const {
  auto out = BinOrEmpty(bin_size, num_bins);
  CHECK(!out.empty());
  return out;
}

std::vector<PcapDataTraceBin> PcapDataTrace::BinOrEmpty(
    std::chrono::milliseconds bin_size, size_t num_bins) const {
  size_t bin_size_ms = bin_size.count();
  std::vector<PcapDataTraceBin> out;
  if (trace_pb_.binned_data_file().empty()) {
    return {};
  }

  AllBinsOrDie([bin_size_ms, num_bins, &out](const PBBinnedData& binned_data) {
    CHECK(binned_data.bytes_binned_size() ==
          binned_data.expected_flows_binned_size());
    CHECK(binned_data.expected_flows_binned_size() ==
          binned_data.packets_binned_size());

    size_t num_bins_masked = num_bins;
    if (num_bins == 0) {
      num_bins_masked = binned_data.bytes_binned_size();
    }

    if (binned_data.bin_size_ms() != bin_size_ms ||
        static_cast<size_t>(binned_data.bytes_binned_size()) <
            num_bins_masked) {
      return;
    }

    if (!out.empty()) {
      return;
    }

    out.resize(num_bins_masked);
    for (size_t i = 0; i < num_bins_masked; ++i) {
      out[i].bytes = binned_data.bytes_binned(i);
      out[i].packets = binned_data.packets_binned(i);
      out[i].flows = binned_data.expected_flows_binned(i);
    }
  });

  return out;
}

std::vector<PcapDataTraceBin> PcapDataTrace::BinOrDie(
    std::chrono::milliseconds bin_size) const {
  return BinOrDie(bin_size, 0);
}

void PcapDataTrace::SinglePassCacheFlows() {
  if (!trace_pb_.tcp_syn_flow_summary_file().empty()) {
    LOG(INFO) << "Already cached";
    return;
  }

  // Need to pick a directory for the file that we will create, will store it
  // in the same directory as the first of the trace files.
  const std::string& first_trace_file = trace_pb_.pcap_files(0);
  std::string dir = ExtractDirName(first_trace_file);
  std::string flow_file_name =
      nc::StrCat(dir, "/flows_summary_", id_.ToString());

  nc::SimTimeEventQueue event_queue;
  std::vector<std::unique_ptr<PcapTraceBinner>> binners;
  std::vector<nc::htsim::PacketObserver*> observers_raw;

  TCPFlowStateRecorder recorder(flow_file_name, kCacheSize, &event_queue);
  observers_raw.emplace_back(&recorder);

  auto packet_gen = GetPacketGenerator(&event_queue);
  std::vector<std::unique_ptr<nc::htsim::BulkPacketSource>> sources;
  sources.emplace_back(std::move(packet_gen));

  ObserverPack observer_pack(observers_raw);
  nc::htsim::BulkPacketGenerator bulk_packet_gen("GenId", std::move(sources),
                                                 &observer_pack, &event_queue);
  bulk_packet_gen.StopQueueWhenDone();
  event_queue.RunAndStopIn(std::chrono::hours(999));

  size_t num_flows = recorder.FinishUp();
  LOG(INFO) << "Parsed " << num_flows << " flows";

  trace_pb_.set_tcp_syn_flow_summary_file(flow_file_name);
  parent_store_->Serialize();
}

void PcapDataTrace::SinglePassCacheBins() {
  std::string current_file_name = trace_pb_.binned_data_file();
  if (!current_file_name.empty()) {
    if (nc::File::Exists(current_file_name)) {
      LOG(INFO) << "Already cached";
      return;
    }
    LOG(ERROR) << "File not found " << current_file_name;
  }

  // Need to pick a directory for the file that we will create, will store it
  // in the same directory as the first of the trace files.
  const std::string& first_trace_file = trace_pb_.pcap_files(0);
  std::string dir = ExtractDirName(first_trace_file);
  std::string bins_file_name =
      nc::StrCat(dir, "/trace_binned_", id_.ToString());

  nc::SimTimeEventQueue event_queue;
  std::vector<std::unique_ptr<PcapTraceBinner>> binners;
  std::vector<nc::htsim::PacketObserver*> observers_raw;

  std::vector<std::chrono::milliseconds> bin_sizes = GetDefaultBinSizes();
  for (std::chrono::milliseconds bin_size : bin_sizes) {
    auto binner = nc::make_unique<PcapTraceBinner>(bin_size, &event_queue);
    observers_raw.emplace_back(binner.get());
    binners.emplace_back(std::move(binner));
  }

  auto packet_gen = GetPacketGenerator(&event_queue);
  std::vector<std::unique_ptr<nc::htsim::BulkPacketSource>> sources;
  sources.emplace_back(std::move(packet_gen));

  ObserverPack observer_pack(observers_raw);
  nc::htsim::BulkPacketGenerator bulk_packet_gen("GenId", std::move(sources),
                                                 &observer_pack, &event_queue);
  bulk_packet_gen.StopQueueWhenDone();
  event_queue.RunAndStopIn(std::chrono::hours(999));

  auto bins_output_stream = GetOutputStream(bins_file_name, false);
  for (const auto& binner : binners) {
    if (binner->bins().empty()) {
      continue;
    }

    PBBinnedData binned_data_pb;
    binned_data_pb.set_bin_size_ms(binner->bin_size().count());
    for (const PcapDataTraceBin& bin : binner->bins()) {
      binned_data_pb.add_bytes_binned(bin.bytes);
      binned_data_pb.add_expected_flows_binned(bin.flows);
      binned_data_pb.add_packets_binned(bin.packets);
    }
    CHECK(WriteDelimitedTo(binned_data_pb, bins_output_stream.get()));
  }

  trace_pb_.set_binned_data_file(bins_file_name);
  parent_store_->Serialize();
}

void PcapDataTrace::BinCache(const std::vector<std::chrono::milliseconds>& bins,
                             std::chrono::milliseconds duration) {
  const std::string& first_trace_file = trace_pb_.pcap_files(0);
  std::string dir = ExtractDirName(first_trace_file);
  std::string file_name = nc::StrCat(dir, "/trace_binned_", id_.ToString());

  nc::SimTimeEventQueue event_queue;
  std::vector<std::unique_ptr<PcapTraceBinner>> binners;
  std::vector<nc::htsim::PacketObserver*> observers_raw;
  for (std::chrono::milliseconds bin_size : bins) {
    size_t bin_count = duration.count() / bin_size.count();

    if (!BinOrEmpty(bin_size, bin_count).empty()) {
      continue;
    }

    LOG(INFO) << "Added binner " << bin_size.count() << "ms";
    auto binner = nc::make_unique<PcapTraceBinner>(bin_size, &event_queue);
    observers_raw.emplace_back(binner.get());
    binners.emplace_back(std::move(binner));
  }

  if (binners.empty()) {
    return;
  }

  auto packet_gen = GetPacketGenerator(&event_queue);
  std::vector<std::unique_ptr<nc::htsim::BulkPacketSource>> sources;
  sources.emplace_back(std::move(packet_gen));

  if (duration.count() == 0) {
    LOG(INFO) << "Will run through the entire trace";
    duration = std::chrono::hours(999);
  }

  ObserverPack observer_pack(observers_raw);
  nc::htsim::BulkPacketGenerator bulk_packet_gen("GenId", std::move(sources),
                                                 &observer_pack, &event_queue);
  bulk_packet_gen.StopQueueWhenDone();
  event_queue.RunAndStopIn(duration + std::chrono::seconds(1));

  auto bins_output_stream = GetOutputStream(file_name, true);
  for (const auto& binner : binners) {
    CHECK(!binner->bins().empty());
    PBBinnedData binned_data_pb;
    binned_data_pb.set_bin_size_ms(binner->bin_size().count());
    for (const PcapDataTraceBin& bin : binner->bins()) {
      binned_data_pb.add_bytes_binned(bin.bytes);
      binned_data_pb.add_expected_flows_binned(bin.flows);
      binned_data_pb.add_packets_binned(bin.packets);
    }
    CHECK(WriteDelimitedTo(binned_data_pb, bins_output_stream.get()));
  }
  bins_output_stream->Close();

  trace_pb_.set_binned_data_file(file_name);
  parent_store_->Serialize();
}

double PcapDataTrace::FractionOfVolumeLocal(
    std::chrono::milliseconds threshold) {
  for (const PBDataTraceLocalFraction& local_fraction_pb :
       trace_pb_.local_fractions()) {
    if (local_fraction_pb.threshold_ms() == threshold.count()) {
      return local_fraction_pb.fraction_local();
    }
  }

  double local_total = 0.0;
  double total = 0.0;
  TCPSYNFlowsOrDie(
      [&local_total, &total, threshold](const TCPSYNFlowSummary& flow) {
        uint64_t rtt_nanos = flow.syn_rtt_estimate_nanos();
        uint64_t bytes = flow.bytes_seen();
        double rtt_ms = rtt_nanos / 1000.0 / 1000.0;
        if (rtt_ms < threshold.count()) {
          local_total += bytes;
        }
        total += bytes;
      });
  double fraction_local = local_total / total;

  PBDataTraceLocalFraction* new_local_fraction_pb =
      trace_pb_.add_local_fractions();
  new_local_fraction_pb->set_fraction_local(fraction_local);
  new_local_fraction_pb->set_threshold_ms(threshold.count());
  parent_store_->Serialize();
  return fraction_local;
}

std::vector<std::string> PcapDataTrace::trace_files() const {
  std::vector<std::string> return_vector;
  for (const std::string& file : trace_pb_.pcap_files()) {
    return_vector.emplace_back(file);
  }
  return return_vector;
}

void PcapDataTrace::set_tag(const std::string& tag) {
  trace_pb_.set_tag(tag);
  parent_store_->Serialize();
}

bool operator<(const TraceSubdivision& a, const TraceSubdivision& b) {
  return std::tie(a.n_, a.active_indices_) < std::tie(b.n_, b.active_indices_);
}

std::string TraceSubdivision::ToString() const {
  std::vector<std::string> indices_str;
  for (size_t i = 0; i < active_indices_.size(); ++i) {
    if (active_indices_[i]) {
      indices_str.emplace_back(std::to_string(i));
    }
  }

  return nc::Substitute("$0_$1", n_, nc::Join(indices_str, "_"));
}

PBTraceSubdivision TraceSubdivision::ToProtobuf() const {
  PBTraceSubdivision out;
  out.set_downscale_n(n_);
  for (size_t i = 0; i < active_indices_.size(); ++i) {
    if (active_indices_[i]) {
      out.add_downscale_indices(i);
    }
  }

  return out;
}

TraceSubdivision::TraceSubdivision(size_t n, std::vector<bool> active_indices)
    : n_(n), active_indices_(active_indices) {
  CHECK(n_ > 1);
  CHECK(active_indices_.size() == n_);
  CHECK(std::any_of(active_indices_.begin(), active_indices_.end(),
                    [](bool v) { return v; }));
}

TraceSubdivision::TraceSubdivision(const PBTraceSubdivision& subdivision)
    : n_(subdivision.downscale_n()) {
  active_indices_.resize(n_, false);
  for (size_t index : subdivision.downscale_indices()) {
    active_indices_[index] = true;
  }
}

TraceId::TraceId(const std::string& location, size_t year, size_t month,
                 size_t day, const std::vector<TraceSubdivision>& subdivisions)
    : location_(location),
      year_(year),
      month_(month),
      day_(day),
      subdivisions_(subdivisions) {}

TraceId::TraceId(const PBTraceId& trace_id)
    : location_(trace_id.location()),
      year_(trace_id.year()),
      month_(trace_id.month()),
      day_(trace_id.day()) {
  for (const PBTraceSubdivision& trace_subdivision : trace_id.subdivisions()) {
    subdivisions_.emplace_back(trace_subdivision);
  }
}

std::string TraceId::ToStringHumanReadable() const {
  std::string subdivisions_str;
  std::function<std::string(const TraceSubdivision&)> f = [](
      const TraceSubdivision& trace_subdivision) {
    return trace_subdivision.ToString();
  };
  nc::Join(subdivisions_.begin(), subdivisions_.end(), "_", f,
           &subdivisions_str);

  return nc::Substitute("$0_$1_$2_$3_$4_", location_, month_, day_, year_,
                        subdivisions_str);
}

std::string TraceId::ToString() const {
  nc::MD5 md5(ToStringHumanReadable());
  return md5.hexdigest();
}

std::string TraceId::LocationAndTimestampToString() const {
  return nc::Substitute("$0_$1_$2_$3", location_, month_, day_, year_);
}

PBTraceId TraceId::ToProtobuf() const {
  PBTraceId out;
  out.set_location(location_);
  out.set_year(year_);
  out.set_month(month_);
  out.set_day(day_);
  for (const TraceSubdivision& subdivision : subdivisions_) {
    *out.add_subdivisions() = subdivision.ToProtobuf();
  }

  return out;
}

TraceId TraceId::Subdivide(size_t n,
                           const std::set<size_t>& active_indices) const {
  TraceId new_trace_id = *this;

  std::vector<bool> active_indices_bitmap(n, false);
  for (size_t index : active_indices) {
    active_indices_bitmap[index] = true;
  }

  new_trace_id.subdivisions_.emplace_back(n, active_indices_bitmap);
  return new_trace_id;
}

bool operator<(const TraceId& a, const TraceId& b) {
  return std::tie(a.location_, a.year_, a.month_, a.day_, a.subdivisions_) <
         std::tie(b.location_, b.year_, b.month_, b.day_, b.subdivisions_);
}

std::ostream& operator<<(std::ostream& output, const TraceId& op) {
  output << op.ToString();
  return output;
}

PcapTraceStore::PcapTraceStore(const std::string& file, bool ignore_missing)
    : file_(file) {
  int fd = open(file.c_str(), O_RDONLY);
  if (!ignore_missing) {
    CHECK(fd > 0);
  }

  if (fd > 0) {
    google::protobuf::io::FileInputStream input_stream(fd);
    PBPcapDataTrace trace_pb;
    while (ReadDelimitedFrom(&trace_pb, &input_stream)) {
      TraceId id(trace_pb.id());
      std::unique_ptr<PcapDataTrace>& trace_ptr = traces_[id];
      trace_ptr =
          std::unique_ptr<PcapDataTrace>(new PcapDataTrace(id, trace_pb, this));
      trace_pb.Clear();
    }
    input_stream.Close();
  }
}

std::vector<TraceId> PcapTraceStore::AllIds(
    const std::string& tag_regex) const {
  std::vector<TraceId> out;
  std::regex regex(tag_regex);
  for (const auto& trace_id_and_trace : traces_) {
    const TraceId& trace_id = trace_id_and_trace.first;
    const PcapDataTrace& data_trace = *trace_id_and_trace.second;

    if (std::regex_search(data_trace.tag(), regex)) {
      out.emplace_back(trace_id);
    }
  }

  return out;
}

void PcapTraceStore::Serialize() {
  auto output_stream = GetOutputStream(file_, false);

  for (const auto& trace_id_and_data_trace : traces_) {
    const PcapDataTrace& data_trace = *trace_id_and_data_trace.second;
    WriteDelimitedTo(data_trace.ToProtobuf(), output_stream.get());
  }

  output_stream->Close();
}

static bool AreSameFile(const std::string& path_one,
                        const std::string& path_two) {
  struct stat stat_one;
  struct stat stat_two;

  if (stat(path_one.c_str(), &stat_one) != 0) {
    LOG(ERROR) << "stat failed for " << path_one << ": " << strerror(errno);
    return false;
  }

  if (stat(path_two.c_str(), &stat_two) != 0) {
    LOG(ERROR) << "stat failed for " << path_two << ": " << strerror(errno);
    return false;
  }

  return stat_one.st_dev == stat_two.st_dev &&
         stat_one.st_ino == stat_two.st_ino;
}

static bool AreSameFileList(const std::vector<std::string>& file_list_one,
                            const std::vector<std::string>& file_list_two) {
  if (file_list_one.size() != file_list_two.size()) {
    return false;
  }

  for (size_t i = 0; i < file_list_one.size(); ++i) {
    if (!AreSameFile(file_list_one[i], file_list_two[i])) {
      return false;
    }
  }

  return true;
}

PcapDataTrace* PcapTraceStore::AddTrace(
    const TraceId& id, const std::vector<std::string>& pcap_files) {
  PcapDataTrace* trace = nc::FindSmartPtrOrNull(traces_, id);
  if (trace != nullptr) {
    CHECK(AreSameFileList(trace->trace_files(), pcap_files))
        << "Inconsistent state " << nc::Join(trace->trace_files(), ",")
        << " vs " << nc::Join(pcap_files, ",");
    return trace;
  }

  LOG(INFO) << "Not found trace for " << id << " will construct a new one";

  PBPcapDataTrace new_trace_pb;
  for (const std::string& pcap_file : pcap_files) {
    new_trace_pb.add_pcap_files(pcap_file);
  }

  std::unique_ptr<PcapDataTrace>& trace_ptr = traces_[id];
  trace_ptr =
      std::unique_ptr<PcapDataTrace>(new PcapDataTrace(id, new_trace_pb, this));

  Serialize();
  return trace_ptr.get();
}

PcapDataTrace& PcapTraceStore::GetTraceOrDie(const TraceId& id) {
  return nc::FindSmartPtrOrDie(traces_, id);
}

const PcapDataTrace& PcapTraceStore::GetTraceOrDie(const TraceId& id) const {
  return nc::FindSmartPtrOrDie(traces_, id);
}

std::vector<TraceId> PcapTraceStore::IdsWithMaxRate(
    const std::vector<TraceId>& ids_to_consider, nc::net::Bandwidth rate,
    std::chrono::milliseconds poll_period,
    std::chrono::milliseconds init_window, double multiplier) {
  double rate_remaining_bps = rate.bps();

  std::vector<TraceId> out;
  for (const TraceId id : ids_to_consider) {
    PcapDataTrace& trace = GetTraceOrDie(id);
    AggregateHistory history =
        trace.InitialHistoryOrDie(poll_period, init_window);
    double max_rate_bps = history.max_rate().bps();
    double mean_rate_bps = history.mean_rate().bps();
    double rate_bps =
        mean_rate_bps + (max_rate_bps - mean_rate_bps) * multiplier;

    if (rate_bps < rate_remaining_bps) {
      rate_remaining_bps -= rate_bps;
      out.emplace_back(id);
      continue;
    }

    double fraction = rate_remaining_bps / rate_bps;
    size_t count = fraction * 100;
    if (count == 0) {
      continue;
    }

    std::set<size_t> indices;
    for (size_t i = 0; i < count; ++i) {
      indices.insert(i);
    }

    PcapDataTrace* trace_to_add = nullptr;
    while (true) {
      PcapDataTrace* new_trace = trace.Subdivide(100, indices);
      new_trace->BinCache({poll_period, init_window}, 2 * init_window);
      nc::net::Bandwidth new_trace_max_rate =
          new_trace->InitialHistoryOrDie(poll_period, init_window).max_rate();
      double new_trace_max_rate_bps = new_trace_max_rate.bps();
      if (new_trace_max_rate_bps > rate_remaining_bps) {
        // The new trace's rate overshoots the target, will reduce the bins and
        // try again.
        if (indices.size() == 1) {
          break;
        }

        auto it = indices.end();
        --it;
        indices.erase(it);
        continue;
      }

      trace_to_add = new_trace;
      break;
    }

    if (trace_to_add == nullptr) {
      continue;
    }

    out.emplace_back(trace_to_add->id());
    break;
  }

  return out;
}

}  // namespace e2e
