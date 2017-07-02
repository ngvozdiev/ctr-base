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

static std::chrono::milliseconds GetQueueSizeMs(const nc::htsim::QueueStats& stats, )

void NetInstrument::Record() {
  for (const nc::htsim::Queue* queue : queues_) {
    const nc::htsim::QueueStats& stats = queue->GetStats();


  }
}

}  // namespace ctr
