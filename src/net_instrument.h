#ifndef NET_INSTRUMENT
#define NET_INSTRUMENT

#include <chrono>
#include <vector>

#include "ncode_common/src/htsim/queue.h"
#include "ncode_common/src/htsim/tcp.h"

namespace ctr {

// Instruments the network by periodically recording per-queue and per-link
// stats to metrics.
class NetInstrument : public nc::EventConsumer {
 public:
  static constexpr char kNetInstrumentId[] = "NetInstrument";

  NetInstrument(const std::vector<const nc::htsim::Queue*>& queues,
                const std::vector<nc::htsim::TCPSource*>& tcp_sources,
                std::chrono::milliseconds record_period,
                nc::EventQueue* event_queue);

  void HandleEvent() override;

 private:
  // Records the stats of all queues.
  void Record();

  // How often to poll queue stats.
  const std::chrono::milliseconds period_;

  // The queues to poll.
  std::vector<const nc::htsim::Queue*> queues_;
};

}  // namespace ctr

#endif