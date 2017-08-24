#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/bulk_gen.h"
#include "ncode_common/src/htsim/pcap_consumer.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/viz/grapher.h"
#include "../pcap_data.h"
#include "dist_model.h"
#include "multi_rate_queue.h"

DEFINE_string(pcap_trace_store, "",
              "A file with information about .pcap traces");
DEFINE_uint64(batch_size, 20, "Number of traces in a batch");
DEFINE_uint64(seed, 1, "Seed to use when picking batches");
DEFINE_uint64(num_batches, 20, "Number of batches");
DEFINE_uint64(bin_size_ms, 100, "Bin size in milliseconds");
DEFINE_uint64(duration_ms, 60000, "How long to run each step for");
DEFINE_uint64(num_queues, 20, "How many queues to run");
DEFINE_uint64(num_steps, 1, "How many steps to perform");

namespace ctr {

// A list of bins for a single trace.
using BinList = std::vector<uint64_t>;

// Extracts all bins from a set of traces. Each bin will be equal to
// num_bytes where num_bytes is the number of bytes
// transmitted during the duration of the bin (bin_size).
static std::vector<BinList> GetAllBins(
    const std::vector<PcapDataTrace*>& traces,
    std::chrono::milliseconds bin_size) {
  std::vector<BinList> all_bins;
  for (PcapDataTrace* trace : traces) {
    BinSequence bin_sequence(trace->TracesAndSlices(trace->AllSlices()));
    PcapDataBinCache cache;
    std::vector<TrimmedPcapDataTraceBin> bins =
        bin_sequence.AccumulateBins(bin_size, &cache);

    BinList bins_for_trace;
    bins_for_trace.reserve(bins.size());

    for (const TrimmedPcapDataTraceBin& bin : bins) {
      bins_for_trace.emplace_back(bin.bytes);
    }

    all_bins.emplace_back(bins_for_trace);
  }

  return all_bins;
}

// Gets the bins for a single step from a set of traces. Returns an empty vector
// if a trace is too short.
static std::vector<BinList> GetBinsForStep(const std::vector<BinList>& all_bins,
                                           size_t step_num,
                                           size_t bins_per_step) {
  size_t from = step_num * bins_per_step;
  size_t to = from + bins_per_step;

  std::vector<BinList> out;
  for (const BinList& bin_list : all_bins) {
    if (bin_list.size() < to) {
      return {};
    }

    BinList bins_for_trace;
    bins_for_trace.resize(bins_per_step);
    for (size_t i = from; i < to; ++i) {
      bins_for_trace[i - from] = bin_list[i];
    }

    out.emplace_back(bins_for_trace);
  }

  return out;
}

// N+1 rates evenly distributed, centered +/- spread around the mean rate.
static std::vector<nc::net::Bandwidth> GetRates(
    const std::vector<BinList>& traces, std::chrono::milliseconds duration,
    size_t n, double spread, nc::net::Bandwidth* mean_rate) {
  using namespace nc::net;

  uint64_t total_bytes = 0;
  for (const auto& trace : traces) {
    total_bytes += std::accumulate(trace.begin(), trace.end(), 0ul);
  }

  double seconds = duration.count() / 1000.0;
  *mean_rate = Bandwidth::FromBitsPerSecond(total_bytes * 8 / seconds);
  double mean_rate_mbps = mean_rate->Mbps();

  double min_mbps = mean_rate_mbps - mean_rate_mbps * spread;
  double max_mbps = mean_rate_mbps + mean_rate_mbps * spread;
  double step_mbps = (max_mbps - min_mbps) / n;

  std::vector<Bandwidth> out;
  for (size_t i = 0; i < n; ++i) {
    double rate_mbps = min_mbps + i * step_mbps;
    out.emplace_back(Bandwidth::FromMBitsPerSecond(rate_mbps));
  }

  return out;
}

static void ParseBatch(const std::vector<PcapDataTrace*>& traces,
                       std::chrono::milliseconds bin_size,
                       std::chrono::milliseconds duration, size_t num_queues,
                       size_t batch_num, int max_step_count,
                       nc::viz::NpyArray* info) {
  nc::SimTimeEventQueue event_queue;

  std::vector<std::unique_ptr<nc::htsim::BulkPacketSource>> sources;
  for (const PcapDataTrace* trace : traces) {
    sources.emplace_back(trace->GetPacketGenerator(&event_queue));
  }

  nc::htsim::BulkPacketGenerator generator("PktGen", std::move(sources),
                                           nullptr, &event_queue);

  // Queues are stored here. They are not deallocated immediately becuse there
  // may be events pending for them from the generator.
  std::vector<std::unique_ptr<MultiRateFIFOQueue>> all_queues;

  std::vector<BinList> all_bins = GetAllBins(traces, bin_size);
  size_t bins_per_step = duration.count() / bin_size.count();
  int step_num = -1;
  while (step_num < (max_step_count - 1)) {
    std::vector<BinList> bins_for_step =
        GetBinsForStep(all_bins, ++step_num, bins_per_step);
    if (bins_for_step.empty()) {
      LOG(INFO) << "Done";
      break;
    }

    nc::net::Bandwidth mean_rate;
    std::vector<nc::net::Bandwidth> rates =
        GetRates(bins_for_step, duration, num_queues, 0.2, &mean_rate);
    std::sort(rates.begin(), rates.end());

    auto new_queue = nc::make_unique<MultiRateFIFOQueue>(rates, &event_queue);
    all_queues.emplace_back(std::move(new_queue));
    MultiRateFIFOQueue* multi_queue = all_queues.back().get();
    generator.set_out(multi_queue);

    // Update the rates of the multi-queue and run it for the next 'duration'.
    event_queue.RunAndStopIn(duration);

    std::map<nc::net::Bandwidth, std::unique_ptr<Distribution>> distributions =
        multi_queue->GetDistributions();
    for (const auto& rate_and_distribution : distributions) {
      nc::net::Bandwidth rate = rate_and_distribution.first;
      const Distribution& distribution = *(rate_and_distribution.second);

      std::vector<nc::viz::NpyArray::StringOrNumeric> row = {
          batch_num, traces.size(), step_num, rate.Mbps(), rate == mean_rate};
      for (size_t i = 0; i < 101; ++i) {
        double p = distribution.Percentile(i / 101.0);
        row.emplace_back(p);
      }

      info->AddRow(row);
    }
  }
}

}  // namespace ctr

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  ctr::PcapTraceStore trace_store(FLAGS_pcap_trace_store);
  std::set<ctr::TraceId> all_ids = trace_store.AllIds();

  std::vector<ctr::PcapDataTrace*> traces;
  for (const ctr::TraceId& trace_id : all_ids) {
    ctr::PcapDataTrace& pcap_trace = *trace_store.GetTraceOrNull(trace_id);
    LOG(INFO) << "Parsing " << trace_id;

    traces.emplace_back(&pcap_trace);
  }

  nc::viz::NpyArray::Types types = {
      {"batch_num", nc::viz::NpyArray::UINT32},
      {"traces_in_batch", nc::viz::NpyArray::UINT32},
      {"step_num", nc::viz::NpyArray::UINT32},
      {"rate_mbps", nc::viz::NpyArray::DOUBLE},
      {"is_mean", nc::viz::NpyArray::UINT8}};
  for (uint32_t i = 0; i < 101; ++i) {
    types.emplace_back(nc::StrCat("queue_size_ms_p", i),
                       nc::viz::NpyArray::DOUBLE);
  }
  nc::viz::NpyArray info(types);

  std::mt19937 rnd(FLAGS_seed);
  CHECK(traces.size() >= FLAGS_batch_size) << traces.size() << " vs "
                                           << FLAGS_batch_size;
  for (size_t i = 0; i < FLAGS_num_batches; ++i) {
    std::shuffle(traces.begin(), traces.end(), rnd);

    std::vector<ctr::PcapDataTrace*> batch(FLAGS_batch_size);
    for (size_t j = 0; j < FLAGS_batch_size; ++j) {
      batch[j] = traces[j];
    }

    ctr::ParseBatch(batch, std::chrono::milliseconds(FLAGS_bin_size_ms),
                    std::chrono::milliseconds(FLAGS_duration_ms),
                    FLAGS_num_queues, i, FLAGS_num_steps, &info);
  }

  info.ToDisk("multi_rate_eval_out");
}
