#include "net_instrument.h"

#include <stddef.h>
#include <algorithm>
#include <iostream>

#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_common.h"
#include "controller.h"
#include "metrics/metrics.h"

namespace ctr {

constexpr char NetInstrument::kNetInstrumentId[];
constexpr char InputPacketObserver::kUntaggedPathId[];

static auto* queue_size_pkts =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string>("net_queue_size_pkts",
                                                  "Size of a queue in packets",
                                                  "queue identifier");

static auto* queue_size_bytes =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string>("net_queue_size_bytes",
                                                  "Size of a queue in bytes",
                                                  "queue identifier");

static auto* queue_size_ms =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string>(
            "net_queue_size_ms", "Size of a queue in milliseconds",
            "queue identifier");

static auto* queue_bytes_seen =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string>(
            "net_queue_bytes_seen",
            "Number of bytes a queue has seen. Cumulative.",
            "queue identifier");

static auto* queue_pkts_seen =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string>(
            "net_queue_pkts_seen",
            "Number of pkts a queue has seen. Cumulative.", "queue identifier");

static auto* queue_pkts_dropped =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string>(
            "net_queue_pkts_dropped",
            "Number of pkts a queue has dropped. Cumulative.",
            "queue identifier");

static auto* tcp_source_fast_retx =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, uint32_t, uint32_t, uint32_t, uint32_t>(
            "tcp_src_fast_retx", "Fast retransmissions", "IP source",
            "IP destination", "TCP source port", "TCP destination port");

static auto* tcp_source_retx_timeouts =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, uint32_t, uint32_t, uint32_t, uint32_t>(
            "tcp_src_retx_timeout", "Retransmission timeouts", "IP source",
            "IP destination", "TCP source port", "TCP destination port");

static auto* tcp_source_completion_times =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, uint32_t, uint32_t, uint32_t, uint32_t>(
            "tcp_source_completion_times",
            "Completion time for each TCP flow as measured by the source",
            "IP source", "IP destination", "TCP source port",
            "TCP destination port");

static auto* tcp_source_close_count =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, uint32_t, uint32_t, uint32_t, uint32_t>(
            "tcp_source_close_count",
            "Number of times the TCP connection has gone through slow-start",
            "IP source", "IP destination", "TCP source port",
            "TCP destination port");

static auto* kPathBytesMetric =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string>(
            "path_bytes", "Bytes through a path", "Path");

static auto* kPathPktsMetric =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string>(
            "path_pkts", "Packets through a path", "Path");

static std::chrono::milliseconds GetQueueSizeMs(
    const nc::htsim::QueueStats& stats, nc::net::Bandwidth rate) {
  double bytes_per_sec = rate.bps() / 8;
  double size_sec = stats.queue_size_bytes / std::max(1.0, bytes_per_sec);
  return std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::duration<double>(size_sec));
}

void NetInstrument::Record() {
  for (const nc::htsim::Queue* queue : queues_) {
    const nc::htsim::QueueStats& stats = queue->GetStats();
    const std::string& id = queue->id();
    queue_size_bytes->GetHandle(id)->AddValue(stats.queue_size_bytes);
    queue_size_pkts->GetHandle(id)->AddValue(stats.queue_size_pkts);
    queue_bytes_seen->GetHandle(id)->AddValue(stats.bytes_seen);
    queue_pkts_seen->GetHandle(id)->AddValue(stats.pkts_seen);
    queue_pkts_dropped->GetHandle(id)->AddValue(stats.pkts_dropped);

    std::chrono::milliseconds size_ms = GetQueueSizeMs(stats, queue->GetRate());
    queue_size_ms->GetHandle(id)->AddValue(size_ms.count());
  }
}

void NetInstrument::HandleEvent() {
  Record();
  EnqueueIn(event_queue()->ToTime(period_));
}

NetInstrument::NetInstrument(
    const std::vector<const nc::htsim::Queue*>& queues,
    const std::vector<nc::htsim::TCPSource*>& tcp_sources,
    std::chrono::milliseconds record_period, nc::EventQueue* event_queue)
    : nc::EventConsumer(kNetInstrumentId, event_queue),
      period_(record_period),
      queues_(queues) {
  for (nc::htsim::TCPSource* tcp_source : tcp_sources) {
    tcp_source->set_complection_times_callback([tcp_source, event_queue](
        nc::EventQueueTime duration, uint64_t close_count) {
      const nc::net::FiveTuple& five_tuple = tcp_source->five_tuple();
      auto* ct_handle = tcp_source_completion_times->GetHandle(
          five_tuple.ip_src().Raw(), five_tuple.ip_dst().Raw(),
          five_tuple.src_port().Raw(), five_tuple.dst_port().Raw());

      size_t time_ms = event_queue->TimeToRawMillis(duration);
      ct_handle->AddValue(time_ms);

      auto* cc_handle = tcp_source_close_count->GetHandle(
          five_tuple.ip_src().Raw(), five_tuple.ip_dst().Raw(),
          five_tuple.src_port().Raw(), five_tuple.dst_port().Raw());
      cc_handle->AddValue(close_count);
    });

    tcp_source->set_fast_retx_callback([tcp_source](uint64_t seq_num) {
      const nc::net::FiveTuple& five_tuple = tcp_source->five_tuple();
      auto* handle = tcp_source_fast_retx->GetHandle(
          five_tuple.ip_src().Raw(), five_tuple.ip_dst().Raw(),
          five_tuple.src_port().Raw(), five_tuple.dst_port().Raw());
      handle->AddValue(seq_num);
    });

    tcp_source->set_retx_timeout_callback([tcp_source](uint64_t seq_num) {
      const nc::net::FiveTuple& five_tuple = tcp_source->five_tuple();
      auto* handle = tcp_source_retx_timeouts->GetHandle(
          five_tuple.ip_src().Raw(), five_tuple.ip_dst().Raw(),
          five_tuple.src_port().Raw(), five_tuple.dst_port().Raw());
      handle->AddValue(seq_num);
    });
  }

  EnqueueIn(event_queue->ToTime(period_));
}

InputPacketObserver::InputPacketObserver(
    const controller::Controller* controller,
    std::chrono::milliseconds record_period, nc::EventQueue* event_queue)
    : EventConsumer("InputPacketObserver", event_queue),
      period_(record_period),
      controller_(controller) {
  EnqueueIn(event_queue->ToTime(period_));
}

void InputPacketObserver::HandleEvent() {
  Record();
  EnqueueIn(event_queue()->ToTime(period_));
}

void InputPacketObserver::ObservePacket(const nc::htsim::Packet& pkt) {
  CHECK(pkt.tag() != nc::htsim::kWildPacketTag);

  BytesAndPackets& bytes_and_packets = tag_to_per_path_state_[pkt.tag()];
  bytes_and_packets.first += pkt.size_bytes();
  ++bytes_and_packets.second;
}

void InputPacketObserver::Record() {
  for (const auto& tag_and_path_state : tag_to_per_path_state_) {
    nc::htsim::PacketTag tag = tag_and_path_state.first;
    const BytesAndPackets& bytes_and_packets = tag_and_path_state.second;

    const nc::net::Walk* path = controller_->PathForTagOrNull(tag);
    std::string path_str = path == nullptr
                               ? kUntaggedPathId
                               : path->ToStringNoPorts(*controller_->graph());

    kPathBytesMetric->GetHandle(path_str)->AddValue(bytes_and_packets.first);
    kPathPktsMetric->GetHandle(path_str)->AddValue(bytes_and_packets.second);
  }
}

}  // namespace ctr
