#include "pcap_data.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/repeated_field.h>
#include <sys/errno.h>
#include <sys/stat.h>
#include <sys/fcntl.h>
#include <unistd.h>
#include <cstdbool>
#include <cstring>
#include <initializer_list>
#include <numeric>
#include <string>
#include <type_traits>
#include <fstream>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/htsim/bulk_gen.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lru_cache.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/pcap.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/substitute.h"
#include "common.h"

namespace ctr {

constexpr size_t PcapDataTrace::kFlowCacheSize;
constexpr size_t PcapDataBinCache::kCacheLineBinCount;

void TrimmedPcapDataTraceBin::Combine(const TrimmedPcapDataTraceBin& other) {
  uint64_t bytes_prev = bytes;
  uint32_t flows_prev = flows_enter;
  bytes += other.bytes;
  flows_enter += other.flows_enter;
  CHECK(flows_enter >= flows_prev);
  CHECK(bytes >= bytes_prev);
}

void TrimmedPcapDataTraceBin::CombineWithFraction(
    const TrimmedPcapDataTraceBin& other, double fraction) {
  uint64_t bytes_prev = bytes;
  uint32_t flows_prev = flows_enter;
  bytes += other.bytes * fraction;
  flows_enter += other.flows_enter * fraction;
  CHECK(flows_enter >= flows_prev);
  CHECK(bytes >= bytes_prev);
}

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
    const T& entry, google::protobuf::io::FileOutputStream* output_stream,
    size_t* bytes_written = nullptr) {
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

  if (bytes_written != nullptr) {
    *bytes_written +=
        (size + ::google::protobuf::io::CodedOutputStream::VarintSize32(size));
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

#if defined(__APPLE__) && defined(__MACH__)
#define lseek64 lseek
#define open64 open
#endif

static std::unique_ptr<google::protobuf::io::FileInputStream> GetInputStream(
    const std::string& filename, size_t offset) {
  int fd = open(filename.c_str(), O_RDONLY);
  CHECK(fd > 0) << "Bad input file " << filename << ": " << strerror(errno);
  CHECK(lseek64(fd, offset, SEEK_SET) != -1);

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

  TCPFlowStateCache(
      std::unique_ptr<google::protobuf::io::FileOutputStream> file_output,
      size_t cache_size)
      : nc::LRUCache<nc::net::FiveTuple, TCPFlowState,
                     nc::net::FiveTupleHasher>(cache_size),
        count_(0),
        file_output_(std::move(file_output)) {}

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
  TCPFlowStateRecorder(
      std::unique_ptr<google::protobuf::io::FileOutputStream> file_output,
      size_t cache_size, nc::EventQueue* event_queue)
      : total_bytes_(0),
        total_pkts_(0),
        cache_(std::move(file_output), cache_size),
        event_queue_(event_queue) {}

  void ObservePacket(const nc::htsim::Packet& pkt) override {
    total_bytes_ += pkt.size_bytes();
    ++total_pkts_;

    if (pkt.is_tcp()) {
      const nc::htsim::TCPPacket* tcp_packet =
          static_cast<const nc::htsim::TCPPacket*>(&pkt);

      bool syn = tcp_packet->syn();
      bool ack = tcp_packet->ack();
      bool fin = tcp_packet->fin();
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

template <typename T>
static bool ReadDelimitedFrom(T* message,
                              google::protobuf::io::CodedInputStream* input,
                              size_t* bytes_read = nullptr) {
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

  if (bytes_read != nullptr) {
    *bytes_read +=
        (size + ::google::protobuf::io::CodedOutputStream::VarintSize32(size));
  }

  return true;
}

BalancingObserver::BalancingObserver(
    const std::vector<nc::htsim::PacketObserver*>& observers)
    : observers_(observers) {
  using namespace nc::htsim;

  nc::htsim::MatchRuleKey key(kWildPacketTag, kWildDevicePortNumber,
                              {nc::net::FiveTuple::kDefaultTuple});
  dispatch_rule_ = nc::make_unique<MatchRule>(key);

  // Will add N actions to the rule, all with the same weight, this will
  // split incoming traffic N ways.
  for (size_t i = 0; i < observers.size(); ++i) {
    auto action = nc::make_unique<MatchRuleAction>(nc::net::DevicePortNumber(i),
                                                   kWildPacketTag, 1);
    dispatch_rule_->AddAction(std::move(action));
  }
}

void BalancingObserver::ObservePacket(const nc::htsim::Packet& pkt) {
  const nc::htsim::MatchRuleAction* action = dispatch_rule_->ChooseOrNull(pkt);
  CHECK(action != nullptr);

  uint32_t index = action->output_port().Raw();
  CHECK(index < observers_.size());
  observers_[index]->ObservePacket(pkt);
}

static std::vector<size_t> SerializeBinnersToFile(
    const std::vector<std::unique_ptr<PcapTraceBinner>>& binners,
    const std::string& file) {
  auto output_stream = GetOutputStream(file, false);

  std::vector<size_t> out;
  size_t total = 0;
  for (const auto& binner : binners) {
    out.emplace_back(total);
    for (const PcapDataTraceBin& bin : binner->bins()) {
      PBBin bin_pb;
      bin_pb.set_byte_count(bin.bytes);
      bin_pb.set_packet_count(bin.packets);
      bin_pb.set_enter_flow_count(bin.flows_enter);
      bin_pb.set_exit_flow_count(bin.flows_exit);
      CHECK(WriteDelimitedTo(bin_pb, output_stream.get(), &total));
    }
  }

  output_stream->Close();
  return out;
}

// Appends 'file' to the end of 'output_file'.
static void Append(const std::string& output_file, const std::string& file) {
  std::ofstream of(output_file, std::ios_base::binary | std::ios_base::app);
  of.seekp(0, std::ios_base::end);

  std::ifstream if_(file, std::ios_base::binary);
  of << if_.rdbuf();
}

static std::unique_ptr<nc::htsim::PcapPacketGen> GetPacketGeneratorStatic(
    const std::vector<std::string>& files, nc::EventQueue* event_queue) {
  std::unique_ptr<nc::pcap::OfflineSourceProvider> offline_provider =
      nc::make_unique<nc::pcap::DefaultOfflineSourceProvider>(files);
  return nc::make_unique<nc::htsim::PcapPacketGen>(std::move(offline_provider),
                                                   event_queue);
}

static constexpr char kFlowsTmpFile[] = "flows.tmp";
static constexpr char kBinnedDataTmpFile[] = "binned_data.tmp";

void PcapDataTrace::Init(const PBPcapDataTrace& trace_pb,
                         const std::string& output_file) {
  nc::SimTimeEventQueue event_queue;

  size_t split_count = trace_pb.split_count();
  std::chrono::microseconds bin_size(trace_pb.bin_size_micros());

  // Will store the flows and the binned data in two separate temporary files.
  // When done parsing will populate the missing fields in 'trace_pb' and
  // serialize it at the back of 'output_file'. Will then also append the
  // contents of the two temporary files to the end of 'output_file'.
  std::string flows_output_file = kFlowsTmpFile;
  std::string binned_data_output_file = kBinnedDataTmpFile;

  auto flows_output_stream = GetOutputStream(flows_output_file, false);
  auto recorder = nc::make_unique<TCPFlowStateRecorder>(
      std::move(flows_output_stream), kFlowCacheSize, &event_queue);

  // There will be as many binners as there are splits, all of them will have
  // the same bin size.
  std::vector<std::unique_ptr<PcapTraceBinner>> binners;
  std::vector<nc::htsim::PacketObserver*> binners_raw;
  CHECK(split_count > 0);
  for (size_t i = 0; i < split_count; ++i) {
    binners.emplace_back(
        nc::make_unique<PcapTraceBinner>(bin_size, &event_queue));
    binners_raw.emplace_back(binners.back().get());
  }
  BalancingObserver balancing_observer(binners_raw);

  std::vector<std::string> pcap_files;
  for (const std::string& file : trace_pb.pcap_files()) {
    pcap_files.emplace_back(file);
  }
  auto packet_gen = GetPacketGeneratorStatic(pcap_files, &event_queue);
  packet_gen->set_max_inter_packet_gap(
      std::chrono::duration_cast<nc::pcap::Timestamp>(
          std::chrono::milliseconds(10)));

  std::vector<std::unique_ptr<nc::htsim::BulkPacketSource>> sources;
  sources.emplace_back(std::move(packet_gen));

  ObserverPack observer_pack({recorder.get(), &balancing_observer});
  nc::htsim::BulkPacketGenerator bulk_packet_gen("GenId", std::move(sources),
                                                 &observer_pack, &event_queue);
  bulk_packet_gen.StopQueueWhenDone();
  event_queue.RunAndStopIn(std::chrono::hours(999));

  // Flush the flows to the temporary file.
  size_t num_flows = recorder->FinishUp();
  recorder.reset();

  // And the bins.
  std::vector<size_t> offsets =
      SerializeBinnersToFile(binners, binned_data_output_file);
  size_t flows_output_size = nc::File::FileSizeOrDie(flows_output_file);
  size_t binned_data_output_size =
      nc::File::FileSizeOrDie(binned_data_output_file);

  // Now we have to update 'trace_pb', append it to 'output_file' and append
  // 'binned_data_file' and 'flows_file' as well.
  PBPcapDataTrace new_trace_pb = trace_pb;
  new_trace_pb.set_tcp_syn_flow_summary_count(num_flows);
  new_trace_pb.set_bin_count(binners[0]->bins().size());
  new_trace_pb.set_flow_summaries_start_offset(binned_data_output_size);
  new_trace_pb.set_flow_summaries_end_offset(binned_data_output_size +
                                             flows_output_size);
  for (size_t offset : offsets) {
    new_trace_pb.add_binned_data_start_offsets(offset);
  }

  auto output_stream = GetOutputStream(output_file, true);
  CHECK(WriteDelimitedTo(new_trace_pb, output_stream.get()));
  output_stream->Close();
  output_stream.reset();

  Append(output_file, binned_data_output_file);
  Append(output_file, flows_output_file);
}

static PBPcapDataTrace GetPcapDataTrace(const std::string& file,
                                        size_t offset) {
  PBPcapDataTrace out;
  auto input_stream = GetInputStream(file, offset);
  google::protobuf::io::CodedInputStream coded_input(input_stream.get());
  ReadDelimitedFrom(&out, &coded_input);

  return out;
}

PcapDataTrace::PcapDataTrace(const std::string& file, size_t offset)
    : store_file_(file),
      offset_in_file_(offset),
      trace_pb_(GetPcapDataTrace(file, offset)),
      trace_pb_size_(trace_pb_.ByteSize()),
      id_(trace_pb_.id()) {}

size_t PcapDataTrace::TotalSizeInFile() const {
  return TraceProtobufSize() + trace_pb_.flow_summaries_end_offset();
}

std::unique_ptr<nc::htsim::PcapPacketGen> PcapDataTrace::GetPacketGenerator(
    nc::EventQueue* event_queue) const {
  std::vector<std::string> files;
  for (const std::string& file : trace_pb_.pcap_files()) {
    files.emplace_back(file);
  }

  return GetPacketGeneratorStatic(files, event_queue);
}

size_t PcapDataTrace::TraceProtobufSize() const {
  uint32_t extra =
      ::google::protobuf::io::CodedOutputStream::VarintSize32(trace_pb_size_);
  return trace_pb_size_ + extra;
}

void PcapDataTrace::TCPSYNFlows(
    std::function<void(const TCPSYNFlowSummary& syn_flow_summay)> callback)
    const {
  if (!callback) {
    return;
  }

  size_t count = trace_pb_.tcp_syn_flow_summary_count();
  size_t offset = offset_in_file_ + TraceProtobufSize() +
                  trace_pb_.flow_summaries_start_offset();
  auto file_input = GetInputStream(store_file_, offset);

  google::protobuf::io::CodedInputStream coded_input(file_input.get());
  TCPSYNFlowSummary flow_summary;
  for (size_t i = 0; i < count; ++i) {
    CHECK(ReadDelimitedFrom(&flow_summary, &coded_input));
    callback(flow_summary);
    flow_summary.Clear();
  }
  file_input->Close();
}

void PcapDataTrace::BinsFromDisk(
    TraceSliceIndex slice, size_t start_bin, size_t end_bin,
    std::function<void(const PBBin& bin)> callback) const {
  if (!callback) {
    return;
  }

  size_t offset = offset_in_file_ + TraceProtobufSize() +
                  trace_pb_.binned_data_start_offsets(slice);
  auto file_input = GetInputStream(store_file_, offset);

  PBBin bin_pb;
  google::protobuf::io::CodedInputStream coded_input(file_input.get());
  for (size_t i = 0; i < end_bin; ++i) {
    CHECK(ReadDelimitedFrom(&bin_pb, &coded_input));
    if (i >= start_bin) {
      callback(bin_pb);
    }
    bin_pb.Clear();
  }
  file_input->Close();
}

std::pair<std::vector<TrimmedPcapDataTraceBin>::const_iterator,
          std::vector<TrimmedPcapDataTraceBin>::const_iterator>
PcapDataBinCache::Bins(const PcapDataTrace* trace, TraceSliceIndex slice,
                       size_t start_bin, size_t end_bin) {
  size_t bin_count = trace->ToProtobuf().bin_count();
  size_t range = end_bin - start_bin;
  CHECK(end_bin <= bin_count);
  CHECK(start_bin <= end_bin);

  CachedTrace& cached_trace = cached_bins_[trace][slice];
  std::vector<TrimmedPcapDataTraceBin>& bins = cached_trace.bins;

  size_t last_cached_bin = cached_trace.starting_bin + bins.size();
  if (cached_trace.starting_bin <= start_bin && last_cached_bin >= end_bin) {
    size_t from = start_bin - cached_trace.starting_bin;
    size_t to = from + range;
    return {bins.begin() + from, bins.begin() + to};
  }

  size_t to_bring_in = std::max(kCacheLineBinCount, range);
  bins.clear();
  bins.reserve(to_bring_in);

  size_t to_bring_in_end = std::min(start_bin + to_bring_in, bin_count);
  trace->BinsFromDisk(slice, start_bin, to_bring_in_end,
                      [&bins](const PBBin& bin) { bins.emplace_back(bin); });

  cached_trace.starting_bin = start_bin;
  return {bins.begin(), bins.begin() + range};
}

std::string PcapDataBinCache::CacheStats() const {
  uint64_t traces = 0;
  uint64_t slices = 0;
  uint64_t bins = 0;

  for (const auto& trace_and_other : cached_bins_) {
    ++traces;
    for (const auto& slice_and_cache : trace_and_other.second) {
      ++slices;
      bins += slice_and_cache.second->bins.size();
    }
  }

  return nc::StrCat("traces: ", traces, " slices: ", slices, " bins: ", bins);
}

std::vector<std::string> PcapDataTrace::trace_files() const {
  std::vector<std::string> return_vector;
  for (const std::string& file : trace_pb_.pcap_files()) {
    return_vector.emplace_back(file);
  }
  return return_vector;
}

TraceId::TraceId(const std::string& location, size_t year, size_t month,
                 size_t day)
    : location_(location), year_(year), month_(month), day_(day) {}

TraceId::TraceId(const PBTraceId& trace_id)
    : location_(trace_id.location()),
      year_(trace_id.year()),
      month_(trace_id.month()),
      day_(trace_id.day()) {}

std::string TraceId::ToString() const {
  return nc::Substitute("$0_$1_$2_$3", location_, month_, day_, year_);
}

PBTraceId TraceId::ToProtobuf() const {
  PBTraceId out;
  out.set_location(location_);
  out.set_year(year_);
  out.set_month(month_);
  out.set_day(day_);

  return out;
}

bool operator<(const TraceId& a, const TraceId& b) {
  return std::tie(a.location_, a.year_, a.month_, a.day_) <
         std::tie(b.location_, b.year_, b.month_, b.day_);
}

std::ostream& operator<<(std::ostream& output, const TraceId& op) {
  output << op.ToString();
  return output;
}

PcapTraceStore::PcapTraceStore(const std::string& file) : file_(file) {
  size_t total_offset = 0;
  while (true) {
    auto input_stream = GetInputStream(file, total_offset);
    input_stream->SetCloseOnDelete(true);

    PBPcapDataTrace trace_pb;
    google::protobuf::io::CodedInputStream coded_input(input_stream.get());
    if (!ReadDelimitedFrom(&trace_pb, &coded_input)) {
      break;
    }

    TraceId id(trace_pb.id());
    auto new_trace = nc::make_unique<PcapDataTrace>(file, total_offset);
    total_offset += new_trace->TotalSizeInFile();

    traces_[id] = std::move(new_trace);
  }
}

BinSequence::BinSequence(const std::vector<TraceAndSlice>& traces)
    : traces_(traces) {
  CHECK(!traces_.empty());
}

BinSequence::BinSequence(std::vector<TraceAndSlice>::const_iterator from,
                         std::vector<TraceAndSlice>::const_iterator to)
    : traces_(from, to) {
  CHECK(!traces_.empty());
}

void BinSequence::Combine(const BinSequence& other) {
  CHECK(other.bin_size() == bin_size());

  traces_.insert(traces_.end(), other.traces_.begin(), other.traces_.end());
}

std::pair<uint64_t, uint64_t> BinSequence::TotalBytesAndPackets() const {
  size_t count = bin_count();

  uint64_t bytes = 0;
  uint64_t packets = 0;
  for (const TraceAndSlice& trace_and_slice : traces_) {
    const PcapDataTrace* trace = trace_and_slice.trace;
    TraceSliceIndex slice = trace_and_slice.slice;
    size_t end_bin = trace_and_slice.start_bin + count;

    trace->BinsFromDisk(
        slice, trace_and_slice.start_bin, end_bin,
        [&bytes, &packets, &trace_and_slice](const PBBin& bin_pb) {
          bytes += bin_pb.byte_count() * trace_and_slice.precise_split;
          packets += bin_pb.packet_count() * trace_and_slice.precise_split;
        });
  }

  return {bytes, packets};
}

nc::net::Bandwidth BinSequence::MaxRate(PcapDataBinCache* cache) const {
  using namespace std::chrono;

  uint32_t max_bytes = 0;
  size_t count = bin_count();
  for (const TraceAndSlice& trace_and_slice : traces_) {
    const PcapDataTrace* trace = trace_and_slice.trace;
    TraceSliceIndex slice = trace_and_slice.slice;
    double fraction = trace_and_slice.precise_split;
    size_t end_bin = trace_and_slice.start_bin + count;

    std::vector<TrimmedPcapDataTraceBin>::const_iterator from;
    std::vector<TrimmedPcapDataTraceBin>::const_iterator to;
    std::tie(from, to) =
        cache->Bins(trace, slice, trace_and_slice.start_bin, end_bin);

    while (from != to) {
      uint32_t bytes_in_bin = from->bytes * fraction;
      max_bytes = std::max(max_bytes, bytes_in_bin);
      ++from;
    }
  }

  std::chrono::microseconds bin_size = this->bin_size();
  double bin_size_seconds = duration_cast<duration<double>>(bin_size).count();
  double bits_per_second = (max_bytes * 8) / bin_size_seconds;
  return nc::net::Bandwidth::FromBitsPerSecond(bits_per_second);
}

PBBinSequence BinSequence::ToProtobuf() const {
  PBBinSequence out;
  for (const TraceAndSlice& trace_and_slice : traces_) {
    PBTraceAndSlice* trace_and_slice_pb = out.add_traces_and_slices();
    *trace_and_slice_pb->mutable_id() =
        trace_and_slice.trace->id().ToProtobuf();
    trace_and_slice_pb->set_slice_index(trace_and_slice.slice);
    trace_and_slice_pb->set_start_bin(trace_and_slice.start_bin);
    trace_and_slice_pb->set_end_bin(trace_and_slice.end_bin);
    trace_and_slice_pb->set_fraction(trace_and_slice.precise_split);
  }

  return out;
}

const std::chrono::microseconds BinSequence::bin_size() const {
  CHECK(!traces_.empty());
  return traces_.front().trace->base_bin_size();
}

std::chrono::milliseconds BinSequence::SimulateQueue(
    nc::net::Bandwidth rate, PcapDataBinCache* cache) const {
  double rate_Bps = rate.bps() / 8.0;
  double bins_in_second =
      1.0 / std::chrono::duration<double>(bin_size()).count();
  double bytes_per_bin = rate_Bps / bins_in_second;

  double max_queue_size_bytes = 0;
  std::vector<TrimmedPcapDataTraceBin> bins = AccumulateBinsPrivate(1, cache);
  double queue_size_bytes = 0;
  for (size_t i = 0; i < bins.size(); ++i) {
    max_queue_size_bytes = std::max(max_queue_size_bytes, queue_size_bytes);

    queue_size_bytes += bins[i].bytes;
    queue_size_bytes -= bytes_per_bin;

    queue_size_bytes = std::max(queue_size_bytes, 0.0);
  }

  double max_queue_size_sec = max_queue_size_bytes / rate_Bps;
  return std::chrono::milliseconds(
      static_cast<uint64_t>(max_queue_size_sec * 1000));
}

AggregateHistory BinSequence::GenerateHistory(
    std::chrono::milliseconds history_bin_size, uint64_t flow_count,
    PcapDataBinCache* cache) const {
  using namespace std::chrono;
  std::vector<TrimmedPcapDataTraceBin> bins_combined =
      AccumulateBins(history_bin_size, cache);

  std::vector<uint64_t> bins_for_history;
  for (const auto& bin_combined : bins_combined) {
    bins_for_history.emplace_back(bin_combined.bytes);
  }

  return {bins_for_history, history_bin_size, flow_count};
}

std::vector<std::unique_ptr<BinSequence>> BinSequence::SplitOrDie(
    const std::vector<double>& fractions) const {
  double total = std::accumulate(fractions.begin(), fractions.end(), 0.0);
  CHECK(total <= 1);
  std::vector<std::unique_ptr<BinSequence>> out;

  double cumulative = 0;
  size_t i = 0;
  for (double fraction : fractions) {
    std::vector<TraceAndSlice> new_traces;
    cumulative += fraction;

    for (; i < cumulative * traces_.size(); ++i) {
      new_traces.emplace_back(traces_[i]);
    }

    if (new_traces.empty()) {
      out.emplace_back();
    } else {
      out.emplace_back(nc::make_unique<BinSequence>(new_traces));
    }
  }

  return out;
}

std::vector<std::unique_ptr<BinSequence>> BinSequence::PreciseSplitOrDie(
    const std::vector<double>& fractions) const {
  double total = std::accumulate(fractions.begin(), fractions.end(), 0.0);
  CHECK(total <= 1);
  std::vector<std::unique_ptr<BinSequence>> out;

  for (double fraction : fractions) {
    std::vector<TraceAndSlice> new_traces = traces_;
    for (auto& trace_and_slice : new_traces) {
      trace_and_slice.precise_split *= fraction;
    }
    CHECK(!new_traces.empty());
    out.emplace_back(nc::make_unique<BinSequence>(new_traces));
  }

  return out;
}

std::unique_ptr<BinSequence> BinSequence::CutFromStart(
    size_t offset_from_start) const {
  size_t count = bin_count();
  CHECK(count >= offset_from_start) << count << " vs " << offset_from_start;

  std::vector<TraceAndSlice> new_traces = traces_;
  for (auto& trace_and_slice : new_traces) {
    trace_and_slice.end_bin = trace_and_slice.start_bin + offset_from_start;
  }

  return nc::make_unique<BinSequence>(new_traces);
}

std::unique_ptr<BinSequence> BinSequence::CutFromStart(
    std::chrono::microseconds duration) const {
  std::chrono::microseconds bin_size_micros = bin_size();
  CHECK(duration.count() % bin_size_micros.count() == 0);
  size_t bin_count = duration.count() / bin_size_micros.count();
  return CutFromStart(bin_count);
}

std::unique_ptr<BinSequence> BinSequence::Offset(size_t offset) const {
  size_t count = bin_count();
  CHECK(count >= offset);

  std::vector<TraceAndSlice> new_traces = traces_;
  for (auto& trace_and_slice : new_traces) {
    trace_and_slice.start_bin += offset;
  }

  return nc::make_unique<BinSequence>(new_traces);
}

std::vector<TrimmedPcapDataTraceBin> BinSequence::AccumulateBins(
    std::chrono::microseconds new_bin_size, PcapDataBinCache* cache) const {
  std::chrono::microseconds base_bin_size = bin_size();

  CHECK(new_bin_size >= base_bin_size);
  size_t base_bins_per_bin = new_bin_size.count() / base_bin_size.count();
  CHECK(new_bin_size.count() % base_bin_size.count() == 0);

  return AccumulateBinsPrivate(base_bins_per_bin, cache);
}

nc::net::Bandwidth BinSequence::MeanRate(PcapDataBinCache* cache) const {
  uint32_t total_bytes = 0;
  size_t count = bin_count();
  for (const TraceAndSlice& trace_and_slice : traces_) {
    const PcapDataTrace* trace = trace_and_slice.trace;
    TraceSliceIndex slice = trace_and_slice.slice;
    double fraction = trace_and_slice.precise_split;
    size_t end_bin = trace_and_slice.start_bin + count;

    std::vector<TrimmedPcapDataTraceBin>::const_iterator from;
    std::vector<TrimmedPcapDataTraceBin>::const_iterator to;
    std::tie(from, to) =
        cache->Bins(trace, slice, trace_and_slice.start_bin, end_bin);

    while (from != to) {
      uint32_t bytes_in_bin = from->bytes;
      total_bytes += bytes_in_bin * fraction;
      ++from;
    }
  }

  double rate_bps =
      total_bytes * 8.0 /
      std::chrono::duration<double>(bin_size() * bin_count()).count();
  return nc::net::Bandwidth::FromBitsPerSecond(rate_bps);
}

std::vector<TrimmedPcapDataTraceBin> BinSequence::AccumulateBinsPrivate(
    size_t bin_size_multiplier, PcapDataBinCache* cache) const {
  size_t count = bin_count();

  std::vector<TrimmedPcapDataTraceBin> out(count);
  for (const TraceAndSlice& trace_and_slice : traces_) {
    const PcapDataTrace* trace = trace_and_slice.trace;
    TraceSliceIndex slice = trace_and_slice.slice;
    size_t start_bin = trace_and_slice.start_bin;
    double fraction = trace_and_slice.precise_split;
    size_t end_bin = trace_and_slice.start_bin + count;

    std::vector<TrimmedPcapDataTraceBin>::const_iterator from;
    std::vector<TrimmedPcapDataTraceBin>::const_iterator to;
    std::tie(from, to) = cache->Bins(trace, slice, start_bin, end_bin);
    CHECK(std::distance(from, to) == static_cast<ssize_t>(count));

    if (fraction == 1) {
      for (size_t i = 0; i < count; ++i) {
        out[i].Combine(*std::next(from, i));
      }
    } else {
      for (size_t i = 0; i < count; ++i) {
        out[i].CombineWithFraction(*std::next(from, i), fraction);
      }
    }
  }

  if (bin_size_multiplier == 1) {
    return out;
  }

  // If the multiplier is >1 will have to combine bins.
  std::vector<TrimmedPcapDataTraceBin> out_combined;
  for (size_t i = 0; i < out.size(); ++i) {
    size_t combined_bin_index = i / bin_size_multiplier;
    out_combined.resize(combined_bin_index + 1);

    TrimmedPcapDataTraceBin& combined_bin = out_combined[combined_bin_index];
    combined_bin.Combine(out[i]);
  }

  return out_combined;
}

const std::chrono::microseconds PcapDataTrace::base_bin_size() const {
  return std::chrono::microseconds(trace_pb_.bin_size_micros());
}

TraceSliceSet PcapDataTrace::AllSlices() const {
  TraceSliceSet out;
  for (size_t i = 0; i < trace_pb_.split_count(); ++i) {
    out.Insert(TraceSliceIndex(i));
  }
  return out;
}

std::vector<BinSequence::TraceAndSlice> PcapDataTrace::TracesAndSlices(
    const TraceSliceSet& slices, std::chrono::milliseconds time_offset) const {
  using namespace std::chrono;
  size_t bin_offset = time_offset.count() /
                      duration_cast<milliseconds>(base_bin_size()).count();
  CHECK(bin_offset < trace_pb_.bin_count());

  std::vector<BinSequence::TraceAndSlice> traces_and_slices;
  for (TraceSliceIndex slice : slices) {
    traces_and_slices.emplace_back(this, slice, bin_offset,
                                   trace_pb_.bin_count(), 1.0);
  }

  return traces_and_slices;
}

std::string PcapDataTrace::Summary() const {
  std::string out;
  nc::SubstituteAndAppend(&out, "Trace of $0 files, first is at $1\n",
                          trace_pb_.pcap_files_size(), trace_pb_.pcap_files(0));
  nc::SubstituteAndAppend(&out,
                          "Bin size $0 microseconds, $1 bins, $2 slices\n",
                          trace_pb_.bin_size_micros(), trace_pb_.bin_count(),
                          trace_pb_.split_count());

  BinSequence bin_sequence(TracesAndSlices(AllSlices()));
  uint64_t total_bytes;
  uint64_t total_packets;
  std::tie(total_bytes, total_packets) = bin_sequence.TotalBytesAndPackets();

  nc::SubstituteAndAppend(&out, "Total $0 bytes, $1 pkts\n", total_bytes,
                          total_packets);
  return out;
}

PcapDataTrace* PcapTraceStore::GetTraceOrNull(const TraceId& id) {
  return nc::FindSmartPtrOrNull(traces_, id);
}

const PcapDataTrace* PcapTraceStore::GetTraceOrNull(const TraceId& id) const {
  return nc::FindSmartPtrOrNull(traces_, id);
}

std::string PcapTraceStore::Summary() const {
  std::vector<std::string> out;
  for (const auto& trace_id_and_trace : traces_) {
    out.emplace_back(trace_id_and_trace.second->Summary());
  }
  return nc::Join(out, "\n");
}

std::unique_ptr<BinSequence> PcapTraceStore::BinSequenceFromProtobufOrDie(
    const PBBinSequence& sequence_pb) const {
  std::vector<BinSequence::TraceAndSlice> traces_and_slices;
  for (const auto& trace_and_slice_pb : sequence_pb.traces_and_slices()) {
    TraceId id(trace_and_slice_pb.id());
    TraceSliceIndex slice(trace_and_slice_pb.slice_index());
    size_t from = trace_and_slice_pb.start_bin();
    size_t to = trace_and_slice_pb.end_bin();
    double fraction = trace_and_slice_pb.fraction();

    const PcapDataTrace* trace_ptr = nc::FindOrDie(traces_, id).get();
    CHECK(trace_ptr->AllSlices().Contains(slice));
    CHECK(from <= to);
    CHECK(to <= trace_ptr->ToProtobuf().bin_count())
        << to << " vs " << trace_ptr->ToProtobuf().bin_count();
    CHECK(fraction <= 1);
    CHECK(fraction > 0);

    traces_and_slices.emplace_back(trace_ptr, slice, from, to, fraction);
  }

  return nc::make_unique<BinSequence>(traces_and_slices);
}

std::unique_ptr<BinSequence> PcapTraceStore::ExtendBinSequence(
    const BinSequence& bin_sequence) const {
  std::vector<BinSequence::TraceAndSlice> extended;
  for (const BinSequence::TraceAndSlice& trace_and_slice :
       bin_sequence.traces()) {
    extended.emplace_back(trace_and_slice);

    uint64_t max_bin_count = trace_and_slice.trace->ToProtobuf().bin_count();
    extended.back().end_bin = max_bin_count;
  }

  return nc::make_unique<BinSequence>(extended);
}

void PcapTraceFitStore::AddToStore(nc::net::Bandwidth rate,
                                   const BinSequence& bin_sequence,
                                   const std::string& output_file,
                                   PcapDataBinCache* cache) {
  std::chrono::milliseconds queue_size =
      bin_sequence.SimulateQueue(rate, cache);
  nc::net::Bandwidth mean_rate = bin_sequence.MeanRate(cache);

  PBTracesToFitRate traces_to_fit_rate_pb;
  *traces_to_fit_rate_pb.mutable_bin_sequence() = bin_sequence.ToProtobuf();
  traces_to_fit_rate_pb.set_max_queue_size_ms(queue_size.count());
  traces_to_fit_rate_pb.set_rate_mbps(rate.Mbps());
  traces_to_fit_rate_pb.set_mean_rate_mbps(mean_rate.Mbps());

  LOG(INFO) << "Adding " << traces_to_fit_rate_pb.DebugString() << " to "
            << output_file;

  auto output_stream = GetOutputStream(output_file, true);
  CHECK(WriteDelimitedTo(traces_to_fit_rate_pb, output_stream.get()));
  output_stream->Close();
  output_stream.reset();
}

PcapTraceFitStore::PcapTraceFitStore(const std::string& file,
                                     const PcapTraceStore* store)
    : store_(store) {
  auto input_stream = GetInputStream(file, 0);
  input_stream->SetCloseOnDelete(true);

  while (true) {
    PBTracesToFitRate traces_to_fit_rate_pb;
    google::protobuf::io::CodedInputStream coded_input(input_stream.get());
    if (!ReadDelimitedFrom(&traces_to_fit_rate_pb, &coded_input)) {
      break;
    }

    nc::net::Bandwidth rate = nc::net::Bandwidth::FromMBitsPerSecond(
        traces_to_fit_rate_pb.rate_mbps());
    LOG(INFO) << "Loaded traces that fit " << rate.Mbps();
    std::unique_ptr<BinSequence> bin_sequence =
        store_->BinSequenceFromProtobufOrDie(
            traces_to_fit_rate_pb.bin_sequence());
    rate_to_bin_sequence_[rate].emplace_back(std::move(bin_sequence));
  }
}

std::unique_ptr<BinSequence> PcapTraceFitStore::GetBinSequence(
    nc::net::Bandwidth bw, nc::net::Bandwidth threshold) {
  CHECK(bw > threshold);
  nc::net::Bandwidth to_look_for = bw - threshold;
  auto it = rate_to_bin_sequence_.upper_bound(to_look_for);
  if (it == rate_to_bin_sequence_.end()) {
    return {};
  }

  // 'it' will be greater or equal to to_look_for.
  nc::net::Bandwidth delta = it->first - to_look_for;
  if (delta > threshold * 2) {
    return {};
  }

  std::unique_ptr<BinSequence> to_return;
  std::vector<std::unique_ptr<BinSequence>>& bins = it->second;
  for (auto bins_it = bins.begin(); bins_it != bins.end(); ++bins_it) {
    std::unique_ptr<BinSequence>& bin_sequence = *bins_it;
    if (bin_sequence) {
      to_return = std::move(bin_sequence);
      bins.erase(bins_it);
      break;
    }
  }

  if (bins.empty()) {
    rate_to_bin_sequence_.erase(it);
  }

  return to_return;
}

static constexpr size_t kExtraPadding = 1000;

// Performs binary search to pick a offset as large as possible into
// traces_and_slices.
static void PickRecursive(
    const std::vector<BinSequence::TraceAndSlice>& traces_and_slices,
    std::chrono::microseconds initial_window, nc::net::Bandwidth target_rate,
    size_t from, size_t to, std::unique_ptr<BinSequence>* out,
    nc::net::Bandwidth* out_rate, PcapDataBinCache* cache) {
  size_t delta = (to - from) / 2;
  if (delta == 0) {
    return;
  }

  size_t offset = from + delta;
  std::unique_ptr<BinSequence> new_sequence;
  if (offset > kExtraPadding) {
    size_t delta = offset - kExtraPadding;
    new_sequence = nc::make_unique<BinSequence>(
        traces_and_slices.begin(), std::next(traces_and_slices.begin(), delta));
  } else {
    double fraction = static_cast<double>(offset) / kExtraPadding;
    auto tmp_sequence = nc::make_unique<BinSequence>(
        traces_and_slices.begin(), std::next(traces_and_slices.begin(), 1));
    new_sequence = std::move(tmp_sequence->PreciseSplitOrDie({fraction})[0]);
  }

  std::unique_ptr<BinSequence> bin_sequence =
      new_sequence->CutFromStart(initial_window);
  std::chrono::milliseconds max_queue_size =
      bin_sequence->SimulateQueue(target_rate, cache);

  if (max_queue_size > std::chrono::milliseconds(10)) {
    //    LOG(ERROR) << "Offset " << offset << " max queue size "
    //               << max_queue_size.count() << "ms";
    PickRecursive(traces_and_slices, initial_window, target_rate, from, offset,
                  out, out_rate, cache);
  } else {
    nc::net::Bandwidth mean_rate = bin_sequence->MeanRate(cache);
    if (*out_rate < mean_rate) {
      *out_rate = mean_rate;
      *out = std::move(bin_sequence);
    }

    //    LOG(ERROR) << "Offset " << offset << " max queue size "
    //               << max_queue_size.count() << "ms mean " << mean_rate.Mbps()
    //               << "Mbps";
    PickRecursive(traces_and_slices, initial_window, target_rate, offset, to,
                  out, out_rate, cache);
  }
}

static std::unique_ptr<BinSequence> PickOrDie(
    const std::vector<BinSequence::TraceAndSlice>& traces_and_slices,
    std::chrono::microseconds initial_window, nc::net::Bandwidth target_rate,
    PcapDataBinCache* cache) {
  nc::net::Bandwidth rate = nc::net::Bandwidth::Zero();
  std::unique_ptr<BinSequence> out;
  PickRecursive(traces_and_slices, initial_window, target_rate, 0,
                kExtraPadding + traces_and_slices.size(), &out, &rate, cache);
  if (!out) {
    // Even the smallest fraction of the first trace was not enough to make the
    // aggregate fit, will use that anyway.
    out = nc::make_unique<BinSequence>(traces_and_slices.begin(),
                                       std::next(traces_and_slices.begin(), 1));
    out = std::move(out->PreciseSplitOrDie({1.0 / kExtraPadding})[0]);
    out = out->CutFromStart(initial_window);
  }

  CHECK(out);

  return out;
}

BinSequenceGenerator::BinSequenceGenerator(
    const std::vector<const PcapDataTrace*>& all_traces,
    const std::vector<std::chrono::milliseconds>& offsets) {
  for (const PcapDataTrace* trace : all_traces) {
    for (std::chrono::milliseconds offset : offsets) {
      std::vector<BinSequence::TraceAndSlice> traces_and_slices =
          trace->TracesAndSlices(trace->AllSlices(), offset);
      all_traces_and_slices_.insert(all_traces_and_slices_.end(),
                                    traces_and_slices.begin(),
                                    traces_and_slices.end());
    }
  }
}

std::unique_ptr<BinSequence> BinSequenceGenerator::Next(
    nc::net::Bandwidth target_rate, std::chrono::microseconds init_window,
    PcapDataBinCache* cache, std::mt19937* rnd) const {
  std::vector<BinSequence::TraceAndSlice> all_traces_and_slices_cpy =
      all_traces_and_slices_;
  std::shuffle(all_traces_and_slices_cpy.begin(),
               all_traces_and_slices_cpy.end(), *rnd);

  return PickOrDie(all_traces_and_slices_cpy, init_window, target_rate, cache);
}

}  // namespace e2e
