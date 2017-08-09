#include <gflags/gflags.h>
#include <algorithm>
#include <chrono>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/pcap.h"
#include "ncode_common/src/viz/grapher.h"
#include "flow_track.h"

DEFINE_string(pcap_file, "", "The file to read from");
DEFINE_string(output, "", "Where to save plot");
DEFINE_uint64(max_flows, 1000000, "Max number of flows to track");

using namespace nc::pcap;
using namespace nc::net;

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  CHECK(!FLAGS_pcap_file.empty());
  CHECK(!FLAGS_output.empty());

  std::unique_ptr<OfflineSourceProvider> offline_source_provider =
      nc::make_unique<DefaultOfflineSourceProvider>(FLAGS_pcap_file);

  nc::FlowTracker tracker(FLAGS_max_flows);
  OfflinePcap offline_pcap(std::move(offline_source_provider), &tracker);
  offline_pcap.Run();

  std::vector<double> fractions;

  std::vector<nc::BidirectionalTCPFlowState> bidirectional_flows =
      tracker.BidirectionalTCPFlows();
  for (const nc::BidirectionalTCPFlowState& flow : bidirectional_flows) {
    Timestamp rtt_estimate = flow.EstimateRTT();
    if (rtt_estimate == Timestamp::max()) {
      continue;
    }

    std::vector<nc::DataCluster> clusters = flow.BreakDown(rtt_estimate * 2);
    for (const auto& cluster : clusters) {
      Timestamp cluster_duration = cluster.last_byte_at - cluster.first_byte_at;
      double fraction =
          cluster_duration.count() / static_cast<double>(rtt_estimate.count());
      fractions.emplace_back(std::max(1.0, fraction));
    }
  }

  nc::viz::DataSeries1D data;
  data.data = std::move(fractions);

  nc::viz::PythonGrapher python_grapher(FLAGS_output);
  python_grapher.PlotCDF({}, {data});
}
