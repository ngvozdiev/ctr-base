#ifndef FLOW_COUNTER_PCAP_DATA_H
#define FLOW_COUNTER_PCAP_DATA_H

#include <map>

#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/net/net_common.h"

namespace ctr {

class FlowCounter {
 public:
  FlowCounter(nc::EventQueue* event_queue)
      : period_start_(event_queue->CurrentTime()), event_queue_(event_queue) {}

  // Adds a new packet. This assumes that now is
  void NewPacket(const nc::net::FiveTuple& five_tuple);

  // Returns an estimate for the expected number of flows that the counter has
  // seen since the last call to EstimateCount. Same as calling
  // EstaimteCountConst() and Clear().
  double EstimateCount();

  // A version of EstimateCount that does not automatically clear current
  // flows. Should call Clear() later.
  double EstimateCountConst() const;

  void Clear() { flows_.clear(); }

 private:
  struct FirstAndLast {
    FirstAndLast(nc::EventQueueTime first) : first(first), last(first) {}

    nc::EventQueueTime first;
    nc::EventQueueTime last;
  };

  // Time the last call to EstimateCount occurred (or object construction).
  nc::EventQueueTime period_start_;

  // Records the timestamps of first and last packets seen from flows.
  std::map<nc::net::FiveTuple, FirstAndLast> flows_;

  // The event queue -- used for getting timestamps.
  nc::EventQueue* event_queue_;
};

}  // namespace ctr

#endif
