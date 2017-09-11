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
#include "ncode_common/src/thread_runner.h"
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
DEFINE_uint64(threads, 4, "Number of threads");
DEFINE_uint64(passes, 100, "Number of passes to perform for each rate");

using Input =
    std::tuple<nc::net::Bandwidth, size_t, const ctr::BinSequenceGenerator*>;

static std::unique_ptr<ctr::BinSequence> HandleMatrixElement(
    const Input& input) {
  nc::net::Bandwidth demand;
  size_t i;
  const ctr::BinSequenceGenerator* sequence_generator;
  std::tie(demand, i, sequence_generator) = input;

  // Each thread get its own cache.
  ctr::PcapDataBinCache bin_cache;

  LOG(INFO) << "Looking for traces to match " << demand.Mbps() << "Mbps";
  std::unique_ptr<ctr::BinSequence> bin_sequence;

  // Will generate FLAGS_passes different bin sequences, and then pick the one
  // with the highest mean.
  nc::net::Bandwidth top_mean = nc::net::Bandwidth::Zero();

  // The RNG.
  std::mt19937 rnd(FLAGS_seed + i);

  for (size_t pass_num = 0; pass_num < FLAGS_passes; ++pass_num) {
    LOG(INFO) << "Pass " << pass_num;
    std::unique_ptr<ctr::BinSequence> bin_sequence_candidate =
        sequence_generator->Next(demand, milliseconds(FLAGS_period_duration_ms),
                                 &rnd);
    if (!bin_sequence_candidate) {
      continue;
    }

    nc::net::Bandwidth mean_rate = bin_sequence_candidate->MeanRate(bin_cache);
    if (mean_rate > top_mean) {
      top_mean = mean_rate;
      bin_sequence = std::move(bin_sequence_candidate);
    }
  }

  CHECK(bin_sequence);
  return bin_sequence;
}

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
  ctr::PcapDataBinCache bin_cache;
  ctr::BinSequenceGenerator bin_sequence_generator(
      trace_store.AllTraces(), {milliseconds(0), seconds(20)}, &bin_cache);

  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
      nc::lp::DemandMatrix::LoadRepetitaOrDie(
          nc::File::ReadFileToStringOrDie(FLAGS_traffic_matrix), node_order,
          &graph);
  demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);

  std::vector<Input> inputs;
  for (size_t i = 0; i < demand_matrix->elements().size(); ++i) {
    nc::net::Bandwidth demand = demand_matrix->elements()[i].demand;
    inputs.emplace_back(demand, i, &bin_cache, &bin_sequence_generator);
  }

  std::vector<std::unique_ptr<ctr::BinSequence>> results =
      nc::RunInParallelWithResult<Input, ctr::BinSequence>(
          inputs, [](const Input& input) { return HandleMatrixElement(input); },
          FLAGS_threads);

  for (size_t i = 0; i < demand_matrix->elements().size(); ++i) {
    nc::net::Bandwidth demand = std::get<0>(inputs[i]);
    ctr::PcapTraceFitStore::AddToStore(demand, *results[i], FLAGS_out_fit_store,
                                       &bin_cache);
  }
}
