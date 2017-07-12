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

class MockNetwork;

class MockDevice : public nc::htsim::DeviceInterface {
 public:
  MockDevice(
      const std::string& id,
      const std::map<nc::htsim::MatchRuleKey, BinSequence>& bin_sequences,
      nc::net::IPAddress ip_address, MockNetwork* mock_network,
      const controller::Controller* controller, nc::EventQueue* event_queue)
      : nc::htsim::DeviceInterface(id, ip_address, event_queue),
        controller_(controller),
        mock_network_(mock_network) {
    for (const auto& key_and_bins : bin_sequences) {
      const nc::htsim::MatchRuleKey& key = key_and_bins.first;
      const BinSequence& bins = key_and_bins.second;
      states_.emplace(key, bins);
    }
  }

  void HandleStateUpdate(const nc::htsim::SSCPAddOrUpdate& update);

  void ReplyToRequest(const nc::htsim::SSCPStatsRequest& request,
                      nc::htsim::SSCPStatsReply* reply);

  void HandlePacket(nc::htsim::PacketPtr pkt) override;

  void HandlePacketFromPort(nc::htsim::Port* input_port,
                            nc::htsim::PacketPtr pkt) override;

 private:
  struct PathState {
    PathState(const nc::net::Walk* path, nc::net::DevicePortNumber output_port,
              nc::htsim::PacketTag tag)
        : path(path), stats(output_port, tag), fraction(0), total_syns(0) {}

    const nc::net::Walk* path;
    nc::htsim::ActionStats stats;
    double fraction;
    size_t total_syns;
  };

  struct AggregateState {
    AggregateState(BinSequence bin_sequence) : bin_sequence(bin_sequence) {}

    // The most recent routes for this aggregate.
    std::map<nc::htsim::PacketTag, PathState> paths;

    // This aggregate's bins.
    BinSequence bin_sequence;
  };

  // Per-aggregate state.
  std::map<nc::htsim::MatchRuleKey, AggregateState> states_;

  // Used to get paths from tags.
  const controller::Controller* controller_;

  // The network.
  MockNetwork* mock_network_;

  friend class MockNetwork;
};

class MockNetwork {
 public:
  MockNetwork(const nc::net::GraphStorage* graph,
              const nc::EventQueue* event_queue)
      : graph_(graph), event_queue_(event_queue) {}

  void AdvanceTime();

  void AddMockDevice(MockDevice* device) { devices_.emplace_back(device); }

 private:
  const nc::net::GraphStorage* graph_;
  const nc::EventQueue* event_queue_;

  std::vector<MockDevice*> devices_;
  nc::EventQueueTime last_advance_time_;
};

class MockDeviceFactory : public controller::DeviceFactory {
 public:
  MockDeviceFactory(const controller::Controller* controller)
      : controller_(controller),
        mock_network_(controller->graph(), controller->event_queue()) {}

  virtual std::unique_ptr<nc::htsim::DeviceInterface> NewDevice(
      const std::string& id, nc::net::IPAddress address,
      nc::EventQueue* event_queue) override {
    auto new_device = nc::make_unique<MockDevice>(
        id, bins_[id], address, &mock_network_, controller_, event_queue);
    mock_network_.AddMockDevice(new_device.get());
    return std::move(new_device);
  }

  void AddBinSequence(const std::string& device_id,
                      const nc::htsim::MatchRuleKey& key,
                      const BinSequence& bin_sequence) {
    std::map<nc::htsim::MatchRuleKey, BinSequence>& bins_for_device =
        bins_[device_id];
    bins_for_device.emplace(key, bin_sequence);
  }

 private:
  // Upon construction each device needs to be initialized with the bin
  // sequences that correspond to aggregates coming out of the device.
  std::map<std::string, std::map<nc::htsim::MatchRuleKey, BinSequence>> bins_;

  const controller::Controller* controller_;

  // A mock network used by MockDevices.
  MockNetwork mock_network_;
};

class MockSimNetwork;

class MockSimDevice : public nc::htsim::Device {
 public:
  MockSimDevice(
      const std::string& id,
      const std::map<nc::htsim::MatchRuleKey, BinSequence>& bin_sequences,
      nc::net::IPAddress ip_address, nc::EventQueue* event_queue)
      : nc::htsim::Device(id, ip_address, event_queue) {
    for (const auto& key_and_bins : bin_sequences) {
      const nc::htsim::MatchRuleKey& key = key_and_bins.first;
      const BinSequence& bins = key_and_bins.second;
      states_.emplace(key, bins);
    }
  }

  void HandlePacket(nc::htsim::PacketPtr pkt) override;

  void HandleStateUpdate(const nc::htsim::SSCPAddOrUpdate& update);

  void PostProcessStats(const nc::htsim::SSCPStatsRequest& request,
                        nc::htsim::SSCPStatsReply* reply) override;

 private:
  struct PathState {
    PathState() : fraction(0), bins_cached_from(0) {}

    double fraction;
    nc::CircularArray<uint64_t, 1 << 10> syns;

    // This path's bins. Used as a cache by MockSimNetwork.
    std::vector<PcapDataTraceBin> bins;
    size_t bins_cached_from;
  };

  struct AggregateState {
    AggregateState(BinSequence bin_sequence)
        : rule(nullptr), bin_sequence(bin_sequence) {}

    // The most recent rule on the device.
    nc::htsim::MatchRule* rule;

    // The most recent routes for this aggregate.
    std::map<nc::htsim::PacketTag, PathState> paths;

    // This aggregate's bins.
    BinSequence bin_sequence;
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

  void PrefetchBins(MockSimDevice::AggregateState* aggregate_state,
                    MockSimDevice::PathState* path_state);

  std::vector<MockSimDevice*> devices_;

  size_t last_bin_count_;

  nc::net::DevicePortNumber default_enter_port_;

  nc::EventQueue* event_queue_;
};

class MockSimDeviceFactory : public controller::DeviceFactory {
 public:
  MockSimDeviceFactory(nc::net::DevicePortNumber default_enter_port,
                       nc::EventQueue* event_queue)
      : mock_network_(default_enter_port, event_queue) {}

  virtual std::unique_ptr<nc::htsim::DeviceInterface> NewDevice(
      const std::string& id, nc::net::IPAddress address,
      nc::EventQueue* event_queue) override {
    auto new_device =
        nc::make_unique<MockSimDevice>(id, bins_[id], address, event_queue);
    new_device->set_die_on_fail_to_match(true);
    mock_network_.AddMockDevice(new_device.get());
    return std::move(new_device);
  }

  void AddBinSequence(const std::string& device_id,
                      const nc::htsim::MatchRuleKey& key,
                      const BinSequence& bin_sequence) {
    std::map<nc::htsim::MatchRuleKey, BinSequence>& bins_for_device =
        bins_[device_id];
    bins_for_device.emplace(key, bin_sequence);
  }

  void Init() { mock_network_.Init(); }

 private:
  // Upon construction each device needs to be initialized with the bin
  // sequences that correspond to aggregates coming out of the device.
  std::map<std::string, std::map<nc::htsim::MatchRuleKey, BinSequence>> bins_;

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
