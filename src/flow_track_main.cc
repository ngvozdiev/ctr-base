#include <gflags/gflags.h>
#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/pcap.h"
#include "flow_track.h"

DEFINE_string(pcap_file, "", "The file to read from");
DEFINE_uint64(max_flows, 1000000, "Max number of flows to track");

using namespace nc::pcap;
using namespace nc::net;

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  CHECK(!FLAGS_pcap_file.empty());

  std::unique_ptr<OfflineSourceProvider> offline_source_provider =
      nc::make_unique<DefaultOfflineSourceProvider>(FLAGS_pcap_file);

  nc::FlowTracker tracker(FLAGS_max_flows);
  OfflinePcap offline_pcap(std::move(offline_source_provider), &tracker);
  offline_pcap.Run();

  std::vector<nc::BidirectionalTCPFlowState> bidirectional_flows =
      tracker.BidirectionalTCPFlows();
  for (const nc::BidirectionalTCPFlowState& flow : bidirectional_flows) {
    LOG(INFO) << flow.ToString();
  }
}
