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
DEFINE_double(thin_fraction, 0.0, "Fraction of aggregates to thin");
DEFINE_double(
    thin_magnitude, 0.1,
    "What fraction of slices to keep for aggregates that are thinned");

using Input = std::tuple<nc::net::Bandwidth, size_t, ctr::PcapDataBinCache*,
                         const ctr::BinSequenceGenerator*>;

static std::unique_ptr<ctr::BinSequence> HandleMatrixElement(
    const Input& input) {
  nc::net::Bandwidth demand;
  size_t i;
  ctr::PcapDataBinCache* bin_cache;
  const ctr::BinSequenceGenerator* sequence_generator;
  std::tie(demand, i, bin_cache, sequence_generator) = input;

  // The RNG.
  std::mt19937 rnd(FLAGS_seed + i);

  std::unique_ptr<ctr::BinSequence> bin_sequence = sequence_generator->Next(
      demand, milliseconds(FLAGS_period_duration_ms), bin_cache, &rnd);

  LOG(INFO) << "Looking for traces to match " << demand.Mbps() << "Mbps got "
            << bin_sequence->MeanRate(bin_cache).Mbps();

  CHECK(bin_sequence);
  return bin_sequence;
}

static std::string StripExtension(const std::string& filename) {
  std::size_t found = filename.find_last_of(".");
  return filename.substr(0, found);
}

static std::unique_ptr<ctr::BinSequence> ThinSequence(
    const ctr::BinSequence& bin_sequence, ctr::PcapDataBinCache* cache) {
  nc::net::Bandwidth target_rate = bin_sequence.MeanRate(cache);
  nc::net::Bandwidth max_rate = bin_sequence.MaxRate(cache);

  auto thinned = bin_sequence.Thin(FLAGS_thin_fraction);
  nc::net::Bandwidth bw = thinned->MeanRate(cache);
  double scale = target_rate / bw;
  auto out = std::move(thinned->PreciseSplitOrDie({scale})[0]);

  nc::net::Bandwidth rate_after_thinning = out->MeanRate(cache);
  nc::net::Bandwidth max_rate_after_thinning = out->MaxRate(cache);

  if (rate_after_thinning == nc::net::Bandwidth::Zero()) {
    return bin_sequence.Duplicate();
  }

  LOG(INFO) << "Will thin aggregate with rate " << target_rate.Mbps() << " max "
            << max_rate.Mbps() << " to rate " << rate_after_thinning.Mbps()
            << " max " << max_rate_after_thinning.Mbps();
  return out;
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
      trace_store.AllTracesExcept(filter),
      {milliseconds(0), milliseconds(20000), milliseconds(40000),
       milliseconds(60000), milliseconds(80000), milliseconds(100000),
       milliseconds(120000), milliseconds(140000)});

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

  // Will now thin aggregates if needed.
  std::mt19937 rnd(FLAGS_seed);
  std::shuffle(results.begin(), results.end(), rnd);
  CHECK(FLAGS_thin_fraction <= 1.0);
  CHECK(FLAGS_thin_fraction >= 0.0);

  size_t thin_limit = FLAGS_thin_fraction * results.size();
  for (size_t i = 0; i < thin_limit; ++i) {
    results[i] = ThinSequence(*results[i], &bin_cache);
  }

  for (size_t i = 0; i < demand_matrix->elements().size(); ++i) {
    ctr::PcapTraceFitStore::AddToStore(*results[i], out, &bin_cache);
  }

  return 0;
}
