#include "flow_counter.h"

#include <algorithm>
#include <cstdint>
#include <utility>

namespace ctr {

void FlowCounter::NewPacket(const nc::net::FiveTuple& five_tuple) {
  nc::EventQueueTime now = event_queue_->CurrentTime();
  auto it_and_rest = flows_.emplace(five_tuple, now);
  auto it = it_and_rest.first;
  FirstAndLast& first_and_last_timestamps = it->second;
  first_and_last_timestamps.last = now;
}

double FlowCounter::EstimateCount() {
  double count = EstimateCountConst();
  Clear();
  return count;
}

double FlowCounter::EstimateCountConst() const {
  nc::EventQueueTime period_len = event_queue_->CurrentTime() - period_start_;
  uint64_t period_len_millis = event_queue_->TimeToRawMillis(period_len);

  nc::EventQueueTime total = nc::EventQueueTime::ZeroTime();
  for (const auto& five_tuple_and_times : flows_) {
    const FirstAndLast& first_and_last_timestamp = five_tuple_and_times.second;
    nc::EventQueueTime delta =
        first_and_last_timestamp.last - first_and_last_timestamp.first;
    total += delta;
  }
  uint64_t total_millis = event_queue_->TimeToRawMillis(total);

  uint64_t flow_count = static_cast<double>(total_millis) / period_len_millis;
  return std::max(static_cast<uint64_t>(1), flow_count);
}

}  // namespace ctr
