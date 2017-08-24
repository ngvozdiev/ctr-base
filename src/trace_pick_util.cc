#include <gflags/gflags.h>
#include <stddef.h>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "pcap_data.h"

using namespace std::chrono;

DEFINE_string(topology, "", "The topology");
DEFINE_string(traffic_matrix, "", "A file with a traffic matrix");
DEFINE_string(pcap_trace_store, "trace_store.pb",
              "A file with information about .pcap traces");
DEFINE_string(out_fit_store, "trace_store_fit.pb",
              "A file with a series of PBTraceToFitRate protobufs");
DEFINE_uint64(period_duration_ms, 60000, "Length of the period");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_uint64(seed, 1, "Seed for the RNG");
DEFINE_uint64(passes, 100, "Number of passes to perform for each rate");

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  CHECK(!FLAGS_topology.empty());
  CHECK(!FLAGS_traffic_matrix.empty());

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology), &node_order);
  builder.RemoveMultipleLinks();
  nc::net::GraphStorage graph(builder);

  ctr::PcapTraceStore trace_store(FLAGS_pcap_trace_store);
  ctr::BinSequenceGenerator bin_sequence_generator(
      trace_store.AllTraces(),
      {milliseconds(0), seconds(10), seconds(20), seconds(30)}, FLAGS_seed);

  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
      nc::lp::DemandMatrix::LoadRepetitaOrDie(
          nc::File::ReadFileToStringOrDie(FLAGS_traffic_matrix), node_order,
          &graph);
  demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);

  size_t i = -1;
  for (const auto& matrix_element : demand_matrix->elements()) {
    nc::net::Bandwidth demand = matrix_element.demand;
    LOG(INFO) << "Looking for traces to match " << matrix_element.demand.Mbps()
              << "Mbps";
    std::unique_ptr<ctr::BinSequence> bin_sequence;

    // Will generate FLAGS_passes different bin sequences, and then pick the one
    // with the highest mean.
    nc::net::Bandwidth top_mean = nc::net::Bandwidth::Zero();

    for (size_t pass_num = 0; pass_num < FLAGS_passes; ++pass_num) {
      LOG(INFO) << "Pass " << pass_num;
      std::unique_ptr<ctr::BinSequence> bin_sequence_candidate =
          bin_sequence_generator.Next(demand,
                                      milliseconds(FLAGS_period_duration_ms));
      if (!bin_sequence_candidate) {
        continue;
      }

      nc::net::Bandwidth mean_rate = bin_sequence_candidate->MeanRate();
      if (mean_rate > top_mean) {
        top_mean = mean_rate;
        bin_sequence = std::move(bin_sequence_candidate);
      }
    }

    CHECK(bin_sequence);
    ctr::PcapTraceFitStore::AddToStore(demand, *bin_sequence,
                                       FLAGS_out_fit_store);

    LOG(INFO) << "Bin sequence " << ++i << " / "
              << demand_matrix->elements().size();
  }
}
