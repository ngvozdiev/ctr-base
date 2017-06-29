#ifndef CTR_NET_MOCK_H
#define CTR_NET_MOCK_H

#include <stddef.h>
#include <map>
#include <memory>
#include <string>

#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/htsim/network.h"
#include "ncode_common/src/htsim/packet.h"
#include "ncode_common/src/net/net_common.h"
#include "controller.h"
#include "pcap_data.h"

namespace ctr {
namespace controller {
class Controller;
} /* namespace controller */
} /* namespace ctr */

namespace ctr {

class MockDevice : public nc::htsim::DeviceInterface {
 public:
  MockDevice(const std::string& id, nc::net::IPAddress ip_address,
             const controller::Controller* controller,
             nc::EventQueue* event_queue)
      : nc::htsim::DeviceInterface(id, ip_address, event_queue),
        controller_(controller) {}

  void HandleStateUpdate(const nc::htsim::SSCPAddOrUpdate& update);

  void AdvanceTime();

  void ReplyToRequest(const nc::htsim::SSCPStatsRequest& request,
                      nc::htsim::SSCPStatsReply* reply);

  void HandlePacket(nc::htsim::PacketPtr pkt) override;

  void HandlePacketFromPort(nc::htsim::Port* input_port,
                            nc::htsim::PacketPtr pkt) override;

  void AddAggregate(const nc::htsim::MatchRuleKey& key,
                    BinSequence bin_sequence) {
    states_.emplace(key, bin_sequence);
  }

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

  nc::EventQueueTime last_advance_time_;
};

class MockDeviceFactory : public controller::DeviceFactory {
 public:
  MockDeviceFactory(const controller::Controller* controller)
      : controller_(controller) {}

  virtual std::unique_ptr<nc::htsim::DeviceInterface> NewDevice(
      const std::string& id, nc::net::IPAddress address,
      nc::EventQueue* event_queue) override {
    return nc::make_unique<MockDevice>(id, address, controller_, event_queue);
  }

 private:
  const controller::Controller* controller_;
};

}  // namespace ctr

#endif
