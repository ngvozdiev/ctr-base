#ifndef CTR_MULTI_RATE_QUEUE_H
#define CTR_MULTI_RATE_QUEUE_H

#include <stddef.h>
#include <cstdint>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/packet.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_common.h"
#include "dist_model.h"

namespace ctr {

class MultiRateFIFOQueue;

// A helper class for MultiRateFIFOQueue below.
class ServiceRateHelper : public nc::EventConsumer {
 public:
  ServiceRateHelper(const std::string& id, nc::net::Bandwidth rate,
                    MultiRateFIFOQueue* parent, std::list<uint16_t>* queue,
                    nc::EventQueue* event_queue, bool pop);

  inline nc::EventQueueTime PacketDrainTime(uint16_t bytes) {
    return nc::EventQueueTime(time_per_bit_ * bytes * 8);
  }

  void HandleEvent() override;

  void PacketAdded(uint16_t pkt_bytes);

  std::unique_ptr<Distribution> GetDistribution() const;

  nc::net::Bandwidth GetRate() const { return rate_; }

  void Kill() { to_kill_ = true; }

 private:
  uint64_t QueueSizeTimeMs() const {
    uint64_t rate_bytes_per_millisecond = rate_.bps() / 8 / 1000;
    return size_bytes_ / rate_bytes_per_millisecond;
  }

  // If true will pop the element, not just increase the iterator.
  bool pop_;

  // If true will ignore events and incoming packets.
  bool to_kill_;

  MultiRateFIFOQueue* parent_;

  // Shared among multiple ServiceRateHelpers. Each "packet" is represented by
  // its size in bytes.
  std::list<uint16_t>* queue_;

  // The rate that this 'subqueue' runs at.
  nc::net::Bandwidth rate_;

  // Number of bytes in this 'subqueue'
  uint64_t size_bytes_;

  // Time to process a single bit.
  double time_per_bit_;

  // The distribution of queue sizes, updated every packet.
  std::vector<uint64_t> queue_size_counts_;

  // Iterator to the next element in the deque_.
  std::list<uint16_t>::iterator it_;
};

// A queue that services multiple rates simultaneously. Useful for testing and
// benchmarking.
class MultiRateFIFOQueue : public nc::htsim::PacketHandler {
 public:
  MultiRateFIFOQueue(const std::vector<nc::net::Bandwidth>& rates,
                     nc::EventQueue* event_queue);

  void HandlePacket(nc::htsim::PacketPtr pkt) override;

  std::map<nc::net::Bandwidth, std::unique_ptr<Distribution>> GetDistributions()
      const;

  void Kill();

  void DeactivateHelper();

 private:
  // Number of still active helpers.
  size_t active_helpers_;

  // If set to true will wait for helpers to become inactive and delete itself.
  bool to_kill_;

  // Shared among all rates.
  std::list<uint16_t> queue_;

  // A per-rate helper.
  std::vector<std::unique_ptr<ServiceRateHelper>> helpers_;
};

}  // namespace ctr

#endif
