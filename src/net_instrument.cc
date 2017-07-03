#include "net_instrument.h"

namespace ctr {

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

NetInstrument::NetInstrument(const std::vector<const nc::htsim::Queue*>& queues,
                             std::chrono::milliseconds record_period,
                             nc::EventQueue* event_queue)
    : nc::EventConsumer(kNetInstrumentId, event_queue),
      period_(record_period),
      queues_(queues) {
  EnqueueIn(event_queue->ToTime(period_));
}

}  // namespace ctr
