#include <gflags/gflags.h>
#include <algorithm>
#include <chrono>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/strutil.h"
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
  std::vector<double> cluster_sizes;

  std::vector<nc::BidirectionalTCPFlowState> bidirectional_flows =
      tracker.BidirectionalTCPFlows();
  for (const nc::BidirectionalTCPFlowState& flow : bidirectional_flows) {
    Timestamp rtt_estimate = flow.EstimateRTT();
    if (rtt_estimate == Timestamp::max()) {
      continue;
    }

    rtt_estimate =
        std::max(rtt_estimate, Timestamp(std::chrono::milliseconds(1)));
    rtt_estimate =
        std::min(rtt_estimate, Timestamp(std::chrono::milliseconds(500)));

    std::vector<nc::DataCluster> clusters = flow.BreakDown(rtt_estimate * 5);
    for (const auto& cluster : clusters) {
      Timestamp cluster_duration = cluster.last_byte_at - cluster.first_byte_at;
      bool all_no_payload = (cluster.bytes == cluster.packets * 40);
      if (all_no_payload) {
        continue;
      }

      double fraction =
          cluster_duration.count() / static_cast<double>(rtt_estimate.count());
      fraction = std::max(1.0, fraction);
      fractions.emplace_back(std::ceil(fraction));
      cluster_sizes.emplace_back(cluster.bytes);
    }
  }

  nc::viz::CDFPlot fractions_plot;
  fractions_plot.AddData("", &fractions);
  fractions_plot.PlotToDir(nc::StrCat(FLAGS_output, "_rtt_fractions"));

  nc::viz::CDFPlot cluster_sizes_plot;
  cluster_sizes_plot.AddData("", &cluster_sizes);
  cluster_sizes_plot.PlotToDir(nc::StrCat(FLAGS_output, "_cluster_sizes"));
}
