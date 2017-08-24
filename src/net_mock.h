#ifndef CTR_NET_MOCK_H
#define CTR_NET_MOCK_H

#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ncode_common/src/circular_array.h"
#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/htsim/network.h"
#include "ncode_common/src/htsim/packet.h"
#include "ncode_common/src/net/net_common.h"
#include "controller.h"
#include "pcap_data.h"

namespace ctr {

class MockSimNetwork;

class MockSimDevice : public nc::htsim::Device {
 public:
  MockSimDevice(const std::string& id,
                std::map<nc::htsim::MatchRuleKey,
                         std::unique_ptr<BinSequence>>&& bin_sequences,
                nc::net::IPAddress ip_address, nc::EventQueue* event_queue)
      : nc::htsim::Device(id, ip_address, event_queue) {
    for (auto& key_and_bins : bin_sequences) {
      const nc::htsim::MatchRuleKey& key = key_and_bins.first;
      std::unique_ptr<BinSequence>& bins = key_and_bins.second;
      states_.emplace(key, std::move(bins));
    }
  }

  void HandlePacket(nc::htsim::PacketPtr pkt) override;

  void HandleStateUpdate(const nc::htsim::SSCPAddOrUpdate& update);

  void PostProcessStats(const nc::htsim::SSCPStatsRequest& request,
                        nc::htsim::SSCPStatsReply* reply) override;

 private:
  struct PathState {
    PathState() : bins_cached_from(0) {}

    // This sequence is a subset of the aggregate's bin sequence (in
    // AggregateState) that should go over this path.
    std::unique_ptr<BinSequence> bin_sequence;

    nc::CircularArray<uint64_t, 1 << 10> syns;

    // This path's bins. Used as a cache by MockSimNetwork.
    std::vector<TrimmedPcapDataTraceBin> bins;
    size_t bins_cached_from;
  };

  struct AggregateState {
    AggregateState(std::unique_ptr<BinSequence> initial_bin_sequence)
        : rule(nullptr),
          initial_bin_sequence(std::move(initial_bin_sequence)) {}

    // The most recent rule on the device.
    nc::htsim::MatchRule* rule;

    // The most recent routes for this aggregate.
    std::map<nc::htsim::PacketTag, PathState> paths;

    // This aggregate's bins.
    std::unique_ptr<BinSequence> initial_bin_sequence;
  };

  // Per-aggregate state.
  std::map<nc::htsim::MatchRuleKey, AggregateState> states_;
  friend class MockSimNetwork;
};

class MockSimNetwork : public nc::EventConsumer {
 public:
  MockSimNetwork(nc::net::DevicePortNumber default_enter_port,
                 nc::EventQueue* event_queue)
      : nc::EventConsumer("MockSimNet", event_queue),
        last_bin_count_(0),
        default_enter_port_(default_enter_port),
        event_queue_(event_queue) {}

  void AddMockDevice(MockSimDevice* device) { devices_.emplace_back(device); }

  void Init() { EnqueueNext(); }

  void HandleEvent() override {
    AdvanceTimeToNextBin();
    EnqueueNext();
  }

  PcapDataBinCache* bin_cache() { return &bin_cache_; }

 private:
  static constexpr size_t kPrefetchSize = 100;

  void EnqueueNext() {
    std::chrono::microseconds bin_size = GetBinSize();
    nc::EventQueueTime bin_size_simtime = event_queue()->ToTime(bin_size);
    EnqueueIn(bin_size_simtime);
  }

  std::chrono::microseconds GetBinSize();

  void AdvanceTimeToNextBin();

  nc::htsim::PacketPtr GetDummyPacket(uint32_t size);

  void PrefetchBins(MockSimDevice::PathState* path_state);

  std::vector<MockSimDevice*> devices_;

  size_t last_bin_count_;

  nc::net::DevicePortNumber default_enter_port_;

  nc::EventQueue* event_queue_;

  PcapDataBinCache bin_cache_;
};

class MockSimDeviceFactory : public controller::DeviceFactory {
 public:
  MockSimDeviceFactory(nc::net::DevicePortNumber default_enter_port,
                       nc::EventQueue* event_queue)
      : mock_network_(default_enter_port, event_queue) {}

  virtual std::unique_ptr<nc::htsim::DeviceInterface> NewDevice(
      const std::string& id, nc::net::IPAddress address,
      nc::EventQueue* event_queue) override {
    std::map<nc::htsim::MatchRuleKey, std::unique_ptr<BinSequence>>&
        bins_for_id = bins_[id];

    auto new_device = nc::make_unique<MockSimDevice>(id, std::move(bins_for_id),
                                                     address, event_queue);
    new_device->set_die_on_fail_to_match(true);
    mock_network_.AddMockDevice(new_device.get());
    return std::move(new_device);
  }

  void AddBinSequence(const std::string& device_id,
                      const nc::htsim::MatchRuleKey& key,
                      std::unique_ptr<BinSequence> bin_sequence) {
    std::map<nc::htsim::MatchRuleKey, std::unique_ptr<BinSequence>>&
        bins_for_device = bins_[device_id];
    bins_for_device.emplace(key, std::move(bin_sequence));
  }

  void Init() { mock_network_.Init(); }

  PcapDataBinCache* bin_cache() { return mock_network_.bin_cache(); }

 private:
  // Upon construction each device needs to be initialized with the bin
  // sequences that correspond to aggregates coming out of the device.
  std::map<std::string, std::map<nc::htsim::MatchRuleKey,
                                 std::unique_ptr<BinSequence>>> bins_;

  // A mock network used by MockDevices.
  MockSimNetwork mock_network_;
};

class DefaultDeviceFactory : public controller::DeviceFactory {
 public:
  DefaultDeviceFactory() {}

  virtual std::unique_ptr<nc::htsim::DeviceInterface> NewDevice(
      const std::string& id, nc::net::IPAddress address,
      nc::EventQueue* event_queue) override {
    auto new_device =
        nc::make_unique<nc::htsim::Device>(id, address, event_queue);
    new_device->set_die_on_fail_to_match(true);
    return std::move(new_device);
  }
};

// Gets a flow count from a series of syns.
uint64_t GetFlowCountFromSyns(const std::vector<uint64_t>& syns);

}  // namespace ctr

#endif
