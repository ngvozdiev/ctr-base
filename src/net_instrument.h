#ifndef NET_INSTRUMENT
#define NET_INSTRUMENT

#include <chrono>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include "ncode/htsim/packet.h"
#include "ncode/htsim/queue.h"
#include "ncode/htsim/tcp.h"

namespace ctr {
namespace controller {
class Controller;
} /* namespace controller */
} /* namespace ctr */

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

// A simple component that records per-path packets per period and bits per
// period. It observes all packets that enter the network. It assumes that all
// packets that it sees are tagged with a path's tag. It maintains separate
// metrics for packets with preferential dropping and packets without. Also
// updates the flow recorder.
class InputPacketObserver : public nc::EventConsumer,
                            public nc::htsim::PacketObserver {
 public:
  // Will use this as a path string for paths that have no tag.
  static constexpr char kUntaggedPathId[] = "Untagged";

  InputPacketObserver(const controller::Controller* controller,
                      std::chrono::milliseconds record_period,
                      nc::EventQueue* event_queue);

  void ObservePacket(const nc::htsim::Packet& pkt) override;

  void HandleEvent() override;

 private:
  using BytesAndPackets = std::pair<uint64_t, uint64_t>;

  // Records per-path stats.
  void Record();

  // How often to record path stats.
  const std::chrono::milliseconds period_;

  // Relates packet tags to paths.
  const controller::Controller* controller_;

  // Indexed by tag value.
  std::map<nc::htsim::PacketTag, BytesAndPackets> tag_to_per_path_state_;
};

// Monitors queueing.
class OutputPacketObserver : public nc::htsim::PacketObserver {
 public:
  OutputPacketObserver(const nc::EventQueue* event_queue)
      : event_queue_(event_queue) {}

  void ObservePacket(const nc::htsim::Packet& pkt) override {
    nc::EventQueueTime queueing_time = pkt.queueing_time();
    nc::EventQueueTime prop_time = pkt.propagation_time();
    uint64_t queueing_time_ms = event_queue_->TimeToRawMillis(queueing_time);
    uint64_t prop_time_ms = event_queue_->TimeToRawMillis(prop_time);
    queueing_time_dist_ms_.Add(queueing_time_ms, 1);
    propagation_time_dist_ms_.Add(prop_time_ms, 1);
  }

  const nc::DiscreteDistribution<uint64_t>& queueing_time_dist_ms() const {
    return queueing_time_dist_ms_;
  }

  const nc::DiscreteDistribution<uint64_t>& propagation_time_dist_ms() const {
    return propagation_time_dist_ms_;
  }

  // Records the distribution. Should be called at the end of the simulation.
  void RecordDist();

 private:
  // Keeps track of the queueing/propagation times of all packets that leave the
  // network.
  nc::DiscreteDistribution<uint64_t> queueing_time_dist_ms_;
  nc::DiscreteDistribution<uint64_t> propagation_time_dist_ms_;

  // The event queue, converts to / from sim time.
  const nc::EventQueue* event_queue_;
};

}  // namespace ctr

#endif
