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

constexpr size_t PcapDataTrace::kCacheSize;

void PcapDataTraceBin::Combine(const PcapDataTraceBin& other, double fraction) {
  bytes += other.bytes * fraction;
  packets += other.packets * fraction;
  flows_enter += other.flows_enter * fraction;
  flows_exit += other.flows_exit * fraction;
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
  nc::htsim::MatchRuleAction* action =
      dispatch_rule_->ChooseOrNull(pkt.five_tuple());
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
      std::move(flows_output_stream), kCacheSize, &event_queue);

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
  uint32_t size = trace_pb_.ByteSize();
  uint32_t extra =
      ::google::protobuf::io::CodedOutputStream::VarintSize32(size);
  return size + extra;
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

std::pair<size_t, size_t> FindStartBinOffset(
    const std::map<size_t, size_t>& offsets, size_t bin) {
  auto it = offsets.lower_bound(bin);
  if (it == offsets.end()) {
    return {0, 0};
  }

  if (it->first == bin) {
    return *it;
  }

  if (it == offsets.begin()) {
    return {0, 0};
  }

  --it;
  return *it;
}

void PcapDataTrace::Bins(size_t slice, size_t start_bin, size_t end_bin,
                         std::function<void(const PBBin& bin)> callback) {
  if (!callback) {
    return;
  }

  std::map<size_t, size_t>& hints = bin_start_hints_[slice];
  size_t bin_offset, byte_offset;
  std::tie(bin_offset, byte_offset) = FindStartBinOffset(hints, start_bin);
  size_t offset = offset_in_file_ + TraceProtobufSize() +
                  trace_pb_.binned_data_start_offsets(slice) + byte_offset;
  auto file_input = GetInputStream(store_file_, offset);

  PBBin bin_pb;
  size_t bytes_read = byte_offset;

  google::protobuf::io::CodedInputStream coded_input(file_input.get());
  size_t i = bin_offset;
  for (; i < end_bin; ++i) {
    CHECK(ReadDelimitedFrom(&bin_pb, &coded_input, &bytes_read));
    if (i >= start_bin) {
      callback(bin_pb);
    }
    bin_pb.Clear();
  }
  file_input->Close();

  if (i == end_bin) {
    hints[end_bin] = bytes_read;
    if (hints.size() == kHintsSize) {
      hints.erase(hints.begin());
    }
  }
}

std::vector<PcapDataTraceBin> PcapDataTrace::BinsCombined(
    const std::set<size_t>& slices, size_t start_bin, size_t end_bin) {
  CHECK(end_bin <= trace_pb_.bin_count());

  const CachedBins* cached = nc::FindOrNull(bins_cache_, slices);
  if (cached == nullptr || (cached->from > start_bin || cached->to < end_bin)) {
    uint64_t offset = (end_bin - start_bin) * 10;
    size_t cache_end = std::min(trace_pb_.bin_count(), end_bin + offset);
    size_t cache_start = offset > start_bin ? 0 : start_bin - offset;
    cached = AddToBinCache(slices, cache_start, cache_end);
  }

  CHECK(cached->from <= start_bin && cached->to >= end_bin);
  const std::vector<PcapDataTraceBin>& cached_bins = cached->bins;
  return {cached_bins.begin() + (start_bin - cached->from),
          cached_bins.begin() + (end_bin - cached->from)};
}

const PcapDataTrace::CachedBins* PcapDataTrace::AddToBinCache(
    const std::set<size_t>& slices, size_t from, size_t to) {
  CachedBins& out = bins_cache_[slices];
  out.bins.clear();
  out.from = from;
  out.to = to;

  out.bins.resize(to - from);
  for (size_t slice : slices) {
    size_t i = -1;
    Bins(slice, from, to, [&out, &i](const PBBin& bin_pb) {
      out.bins[++i].Combine(PcapDataTraceBin(bin_pb), 1.0);
    });
  }

  size_t bytes = sizeof(PcapDataTraceBin) * (to - from);
  LOG(ERROR) << "Added to cache " << slices.size() << " slices from " << from
             << " to " << to << " size: " << bytes / 1000.0 / 1000.0 << "MB";
  return &out;
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

void BinSequence::Combine(const BinSequence& other) {
  CHECK(other.bin_size() == bin_size());

  traces_.insert(traces_.end(), other.traces_.begin(), other.traces_.end());
}

std::pair<uint64_t, uint64_t> BinSequence::TotalBytesAndPackets() {
  size_t count = bin_count();

  uint64_t bytes = 0;
  uint64_t packets = 0;
  for (const TraceAndSlice& trace_and_slice : traces_) {
    PcapDataTrace* trace = trace_and_slice.trace;
    size_t slice = trace_and_slice.slice;
    size_t end_bin = trace_and_slice.start_bin + count;
    trace->Bins(slice, trace_and_slice.start_bin, end_bin,
                [&bytes, &packets, &trace_and_slice](const PBBin& bin_pb) {
                  bytes += bin_pb.byte_count() * trace_and_slice.fraction;
                  packets += bin_pb.packet_count() * trace_and_slice.fraction;
                });
  }

  return {bytes, packets};
}

const std::chrono::microseconds BinSequence::bin_size() const {
  CHECK(!traces_.empty());
  return traces_.front().trace->base_bin_size();
}

std::vector<double> BinSequence::Residuals(nc::net::Bandwidth rate) {
  double rate_Bps = rate.bps() / 8.0;
  double bins_in_second =
      1.0 / std::chrono::duration<double>(bin_size()).count();
  double bytes_per_bin = rate_Bps / bins_in_second;

  std::vector<PcapDataTraceBin> bins = AccumulateBinsPrivate(1);
  std::vector<double> out(bins.size());
  for (size_t i = 0; i < bins.size(); ++i) {
    out[i] = (bins[i].bytes - bytes_per_bin) / bytes_per_bin;
  }

  return out;
}

AggregateHistory BinSequence::GenerateHistory(
    std::chrono::milliseconds history_bin_size) {
  using namespace std::chrono;
  std::vector<PcapDataTraceBin> bins_combined =
      AccumulateBins(history_bin_size);

  size_t total_flows = 0;
  std::vector<uint64_t> bins_for_history;
  for (const auto& bin_combined : bins_combined) {
    bins_for_history.emplace_back(bin_combined.bytes);
    total_flows += bin_combined.flows_enter;
  }

  return {bins_for_history, history_bin_size,
          std::max(1ul, total_flows / bins_combined.size())};
}

std::vector<BinSequence> BinSequence::SplitOrDie(
    const std::vector<double>& fractions) const {
  double total = std::accumulate(fractions.begin(), fractions.end(), 0.0);
  CHECK(total <= 1);
  std::vector<BinSequence> out;

  double cumulative = 0;
  size_t i = 0;
  for (double fraction : fractions) {
    std::vector<TraceAndSlice> new_traces;
    cumulative += fraction;

    for (; i < cumulative * traces_.size(); ++i) {
      new_traces.emplace_back(traces_[i]);
    }

    out.emplace_back(new_traces);
  }

  return out;
}

std::vector<BinSequence> BinSequence::PreciseSplitOrDie(
    const std::vector<double>& fractions) const {
  std::vector<BinSequence> out;

  for (double fraction : fractions) {
    std::vector<TraceAndSlice> new_traces = traces_;
    for (auto& trace_and_slce : new_traces) {
      trace_and_slce.fraction *= fraction;
    }

    out.emplace_back(new_traces);
  }

  return out;
}

BinSequence BinSequence::CutFromStart(size_t offset_from_start) const {
  size_t count = bin_count();
  CHECK(count >= offset_from_start) << count << " vs " << offset_from_start;

  std::vector<TraceAndSlice> new_traces = traces_;
  for (auto& trace_and_slice : new_traces) {
    trace_and_slice.end_bin = trace_and_slice.start_bin + offset_from_start;
  }

  return {new_traces};
}

BinSequence BinSequence::Offset(size_t offset) const {
  size_t count = bin_count();
  CHECK(count >= offset);

  std::vector<TraceAndSlice> new_traces = traces_;
  for (auto& trace_and_slice : new_traces) {
    trace_and_slice.start_bin += offset;
  }

  return {new_traces};
}

std::vector<PcapDataTraceBin> BinSequence::AccumulateBins(
    std::chrono::microseconds new_bin_size) {
  std::chrono::microseconds base_bin_size = bin_size();

  CHECK(new_bin_size >= base_bin_size);
  size_t base_bins_per_bin = new_bin_size.count() / base_bin_size.count();
  CHECK(new_bin_size.count() % base_bin_size.count() == 0);

  return AccumulateBinsPrivate(base_bins_per_bin);
}

nc::net::Bandwidth BinSequence::MeanRate() {
  uint64_t total_bytes;
  uint64_t total_packets;
  std::tie(total_bytes, total_packets) = TotalBytesAndPackets();

  double rate_bps =
      total_bytes * 8.0 /
      std::chrono::duration<double>(bin_size() * bin_count()).count();
  return nc::net::Bandwidth::FromBitsPerSecond(rate_bps);
}

std::vector<PcapDataTraceBin> BinSequence::AccumulateBinsPrivate(
    size_t bin_size_multiplier) {
  size_t count = bin_count();
  // Will first group slices by <trace, start_bin, end_bin, fraction>.
  std::map<std::tuple<PcapDataTrace*, size_t, size_t, double>, std::set<size_t>>
      slices_grouped;
  for (const TraceAndSlice& trace_and_slice : traces_) {
    PcapDataTrace* trace = trace_and_slice.trace;
    size_t slice = trace_and_slice.slice;
    size_t start_bin = trace_and_slice.start_bin;
    size_t end_bin = trace_and_slice.start_bin + count;
    slices_grouped[std::make_tuple(trace, start_bin, end_bin,
                                   trace_and_slice.fraction)]
        .insert(slice);
  }

  std::vector<PcapDataTraceBin> out(count);
  for (const auto& trace_and_slices : slices_grouped) {
    PcapDataTrace* trace;
    size_t start_bin, end_bin;
    double fraction;
    std::tie(trace, start_bin, end_bin, fraction) = trace_and_slices.first;
    const std::set<size_t>& slices = trace_and_slices.second;

    std::vector<PcapDataTraceBin> bins_combined =
        trace->BinsCombined(slices, start_bin, end_bin);
    CHECK(bins_combined.size() == count);

    for (size_t i = 0; i < count; ++i) {
      out[i].Combine(bins_combined[i], fraction);
    }
  }

  if (bin_size_multiplier == 1) {
    return out;
  }

  // If the multiplier is >1 will have to combine bins.
  std::vector<PcapDataTraceBin> out_combined;
  for (size_t i = 0; i < out.size(); ++i) {
    size_t combined_bin_index = i / bin_size_multiplier;
    out_combined.resize(combined_bin_index + 1);

    PcapDataTraceBin& combined_bin = out_combined[combined_bin_index];
    combined_bin.bytes += out[i].bytes;
    combined_bin.packets += out[i].packets;
    combined_bin.flows_enter += out[i].flows_enter;
    combined_bin.flows_exit += out[i].flows_exit;
  }

  return out_combined;
}

const std::chrono::microseconds PcapDataTrace::base_bin_size() const {
  return std::chrono::microseconds(trace_pb_.bin_size_micros());
}

std::set<size_t> PcapDataTrace::AllSlices() const {
  std::set<size_t> out;
  for (size_t i = 0; i < trace_pb_.split_count(); ++i) {
    out.emplace(i);
  }
  return out;
}

BinSequence PcapDataTrace::ToSequence(const std::set<size_t>& slices) {
  std::vector<BinSequence::TraceAndSlice> traces_and_slices;
  for (size_t slice : slices) {
    CHECK(slice < trace_pb_.split_count());
    traces_and_slices.push_back({this, slice, 0, trace_pb_.bin_count(), 1.0});
  }

  return {traces_and_slices};
}

std::string PcapDataTrace::Summary() {
  std::string out;
  nc::SubstituteAndAppend(&out, "Trace of $0 files, first is at $1\n",
                          trace_pb_.pcap_files_size(), trace_pb_.pcap_files(0));
  nc::SubstituteAndAppend(&out,
                          "Bin size $0 microseconds, $1 bins, $2 slices\n",
                          trace_pb_.bin_size_micros(), trace_pb_.bin_count(),
                          trace_pb_.split_count());

  BinSequence bin_sequence = ToSequence(AllSlices());
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

BinSequence BinsAtRate(nc::net::Bandwidth target_rate,
                       std::chrono::microseconds init_window,
                       BinSequenceGenerator* sequence_generator) {
  nc::net::Bandwidth rate_remaining = target_rate;

  std::unique_ptr<BinSequence> out;
  while (true) {
    BinSequence sequence = sequence_generator->Next();
    std::chrono::microseconds bin_size = sequence.bin_size();
    size_t bin_count = init_window.count() / bin_size.count();
    CHECK(bin_count > 0);

    BinSequence sub_sequence = sequence.CutFromStart(bin_count);
    LOG(ERROR) << "A";
    nc::net::Bandwidth rate = sub_sequence.MeanRate();
    LOG(ERROR) << "B " << rate.Mbps() << " r " << rate_remaining.Mbps();

    if (rate <= rate_remaining) {
      rate_remaining -= rate;
      if (out) {
        out->Combine(sequence);
      } else {
        out = nc::make_unique<BinSequence>(sequence);
      }

      continue;
    }

    double fraction = rate_remaining / rate;
    BinSequence new_sequence = sequence.SplitOrDie({fraction})[0];
    if (out) {
      out->Combine(new_sequence);
    } else {
      out = nc::make_unique<BinSequence>(new_sequence);
    }

    break;
  }

  return *out;
}

BinSequence BinSequenceGenerator::Next() {
  const BinSequence& trace_to_offset = initial_list_[i_ % initial_list_.size()];
  size_t step_count = i_++ / initial_list_.size();
  return trace_to_offset.Offset(step_count * offset_step_);
}

}  // namespace e2e
