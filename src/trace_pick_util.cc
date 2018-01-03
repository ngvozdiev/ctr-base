#include <gflags/gflags.h>
#include <stddef.h>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "ncode/common.h"
#include "ncode/file.h"
#include "ncode/logging.h"
#include "ncode/strutil.h"
#include "ncode/lp/demand_matrix.h"
#include "ncode/net/net_common.h"
#include "ncode/net/net_gen.h"
#include "ncode/thread_runner.h"
#include "pcap_data.h"

using namespace std::chrono;

DEFINE_string(topology, "", "The topology");
DEFINE_string(traffic_matrix, "", "A file with a traffic matrix");
DEFINE_string(pcap_trace_store, "trace_store.pb",
              "A file with information about .pcap traces");
DEFINE_string(out, "",
              "Will output a series of PBTraceToFitRate protobufs, they will "
              "only be compatible with the same trace_store.pb");
DEFINE_uint64(period_duration_ms, 60000, "Length of the period");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_uint64(seed, 1, "Seed for the RNG");
DEFINE_uint64(passes, 10, "Number of passes to perform for each rate");
DEFINE_string(to_exclude,
              "equinix-sanjose_A_12_19_2013,equinix-chicago_A_11_21_2013",
              "Traces to exclude");

using Input = std::tuple<nc::net::Bandwidth, size_t, ctr::PcapDataBinCache*,
                         const ctr::BinSequenceGenerator*>;

static std::unique_ptr<ctr::BinSequence> HandleMatrixElement(
    const Input& input) {
  nc::net::Bandwidth demand;
  size_t i;
  ctr::PcapDataBinCache* bin_cache;
  const ctr::BinSequenceGenerator* sequence_generator;
  std::tie(demand, i, bin_cache, sequence_generator) = input;

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
                                 bin_cache, &rnd);
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

static std::string StripExtension(const std::string& filename) {
  std::size_t found = filename.find_last_of(".");
  return filename.substr(0, found);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  CHECK(!FLAGS_topology.empty());
  CHECK(!FLAGS_traffic_matrix.empty());

  std::string out = FLAGS_out;
  if (out.empty()) {
    out = nc::StrCat(StripExtension(FLAGS_traffic_matrix), "_fit.pb");
  }

  if (nc::File::Exists(out)) {
    LOG(ERROR) << "Output " << out << " already exists";
    return 0;
  }

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology), &node_order);
  builder.RemoveMultipleLinks();
  nc::net::GraphStorage graph(builder);

  std::vector<std::string> to_exclude_v = nc::Split(FLAGS_to_exclude, ",");
  std::set<std::string> to_exclude(to_exclude_v.begin(), to_exclude_v.end());
  auto filter = ctr::PcapTraceStore::FilterFunction(
      [&to_exclude](const ctr::TraceId& trace_id) {
        if (nc::ContainsKey(to_exclude, trace_id.ToString())) {
          LOG(INFO) << "Filtering " << trace_id.ToString();
          return true;
        }

        return false;
      });

  ctr::PcapTraceStore trace_store(FLAGS_pcap_trace_store);
  ctr::BinSequenceGenerator bin_sequence_generator(
      trace_store.AllTracesExcept(filter), {milliseconds(0)});

  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
      nc::lp::DemandMatrix::LoadRepetitaFileOrDie(FLAGS_traffic_matrix,
                                                  node_order, &graph);
  demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);

  // The bin cache.
  ctr::PcapDataBinCache bin_cache;

  std::vector<Input> inputs;
  for (size_t i = 0; i < demand_matrix->elements().size(); ++i) {
    nc::net::Bandwidth demand = demand_matrix->elements()[i].demand;
    inputs.emplace_back(demand, i, &bin_cache, &bin_sequence_generator);
  }

  // It is tricky to get good-performance thread-safe code in pcap_data.cc
  // because of the copies the cache would have to do. So for now we only do one
  // thread.
  std::vector<std::unique_ptr<ctr::BinSequence>> results =
      nc::RunInParallelWithResult<Input, ctr::BinSequence>(
          inputs, [](const Input& input) { return HandleMatrixElement(input); },
          1);

  for (size_t i = 0; i < demand_matrix->elements().size(); ++i) {
    nc::net::Bandwidth demand = std::get<0>(inputs[i]);
    ctr::PcapTraceFitStore::AddToStore(demand, *results[i], out, &bin_cache);
  }

  return 0;
}
