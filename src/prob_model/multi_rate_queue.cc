#include "multi_rate_queue.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>

#include "ncode/common.h"
#include "ncode/strutil.h"

namespace ctr {

ServiceRateHelper::ServiceRateHelper(const std::string& id,
                                     nc::net::Bandwidth rate,
                                     MultiRateFIFOQueue* parent,
                                     std::list<uint16_t>* queue,
                                     nc::EventQueue* event_queue, bool pop)
    : nc::EventConsumer(id, event_queue),
      pop_(pop),
      to_kill_(false),
      parent_(parent),
      queue_(queue),
      rate_(rate),
      size_bytes_(0) {
  time_per_bit_ = pow(10, 12.0) / static_cast<double>(rate_.bps());
  it_ = queue->end();
}

void ServiceRateHelper::PacketAdded(uint16_t pkt_bytes) {
  CHECK(time_per_bit_ != 0);

  if (to_kill_) {
    return;
  }

  uint64_t size_ms = QueueSizeTimeMs();
  size_ms = std::min(size_ms, static_cast<uint64_t>(1000));
  queue_size_counts_.resize(
      std::max(queue_size_counts_.size(), static_cast<size_t>(size_ms + 1)));
  ++queue_size_counts_[size_ms];

  size_bytes_ += pkt_bytes;
  if (outstanding_event_count() == 0) {
    if (it_ == queue_->end()) {
      it_ = std::prev(queue_->end(), 1);
    }

    EnqueueIn(PacketDrainTime(*it_));
  }
}

std::unique_ptr<Distribution> ServiceRateHelper::GetDistribution() const {
  uint64_t total = std::accumulate(queue_size_counts_.begin(),
                                   queue_size_counts_.end(), 0ul);
  std::vector<double> probabilities(queue_size_counts_.size(), 0.0);
  for (size_t i = 0; i < queue_size_counts_.size(); ++i) {
    probabilities[i] = queue_size_counts_[i] / static_cast<double>(total);
  }

  return nc::make_unique<Distribution>(probabilities, 0ul);
}

MultiRateFIFOQueue::MultiRateFIFOQueue(
    const std::vector<nc::net::Bandwidth>& rates, nc::EventQueue* event_queue)
    : active_helpers_(0), to_kill_(false) {
  CHECK(!rates.empty());
  for (uint32_t i = 0; i < rates.size(); ++i) {
    auto new_ptr = nc::make_unique<ServiceRateHelper>(
        nc::StrCat("R_", i), rates[i], this, &queue_, event_queue, i == 0);
    helpers_.emplace_back(std::move(new_ptr));
  }
}

void MultiRateFIFOQueue::HandlePacket(nc::htsim::PacketPtr pkt) {
  uint16_t size_bytes = pkt->size_bytes();
  queue_.emplace_back(size_bytes);
  for (auto& helper : helpers_) {
    helper->PacketAdded(size_bytes);
  }
}

std::map<nc::net::Bandwidth, std::unique_ptr<Distribution>>
MultiRateFIFOQueue::GetDistributions() const {
  std::map<nc::net::Bandwidth, std::unique_ptr<Distribution>> out;
  for (const auto& queue : helpers_) {
    out.emplace(queue->GetRate(), queue->GetDistribution());
  }

  return out;
}

void MultiRateFIFOQueue::Kill() {
  to_kill_ = true;
  for (auto& helper : helpers_) {
    helper->Kill();
    if (helper->outstanding_event_count() > 0) {
      ++active_helpers_;
    }
  }

  if (active_helpers_ == 0) {
    delete this;
  }
}

void MultiRateFIFOQueue::DeactivateHelper() {
  if (!to_kill_) {
    return;
  }

  CHECK(active_helpers_ > 0);
  --active_helpers_;
  if (active_helpers_ == 0) {
    delete this;
  }
}

void ServiceRateHelper::HandleEvent() {
  if (to_kill_) {
    if (outstanding_event_count() == 0) {
      parent_->DeactivateHelper();
    }

    return;
  }

  CHECK(it_ != queue_->end());
  CHECK(time_per_bit_ != 0);

  uint16_t pkt_bytes = *it_;
  CHECK(size_bytes_ >= pkt_bytes);
  size_bytes_ -= pkt_bytes;

  auto it_prev = it_;
  ++it_;
  if (pop_) {
    queue_->erase(it_prev);
  }

  if (it_ != queue_->end()) {
    EnqueueIn(PacketDrainTime(*it_));
  }
}

}  // namespace ctr
