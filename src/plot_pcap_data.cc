#include <gflags/gflags.h>
#include <stddef.h>
#include <chrono>
#include <string>
#include <vector>

#include "ncode_common/src/stats.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/viz/grapher.h"
#include "pcap_data.h"

DEFINE_string(pcap_trace_store, "trace_store.pb",
              "A file with information about .pcap traces");
DEFINE_string(traces, "", "Traces to plot, if empty will plot all traces");

using namespace std::chrono;

// static std::vector<double> Bin(const std::vector<double>& data,
//                               size_t bin_size) {
//  size_t max_count = (data.size() / bin_size) * bin_size;
//
//  std::vector<double> out;
//  for (size_t i = 0; i < max_count; ++i) {
//    size_t bin_index = i / bin_size;
//    out.resize(bin_index + 1);
//    out[bin_index] += data[bin_index];
//  }
//
//  for (double& v : out) {
//    v /= bin_size;
//  }
//
//  return out;
//}
//
// static std::vector<std::pair<double, double>> AllanVariance(
//    const std::vector<double>& data) {
//  std::vector<std::pair<double, double>> out;
//  for (size_t bin_size = 1; bin_size < data.size() / 2; ++bin_size) {
//    LOG(INFO) << "AV bin size " << bin_size;
//    std::vector<double> binned = Bin(data, bin_size);
//
//    nc::SummaryStats stats;
//    for (double bin : binned) {
//      stats.Add(bin);
//    }
//
//    out.emplace_back(bin_size, stats.var());
//  }
//
//  return out;
//}

static std::vector<double> BinTrace(const ctr::BinSequence& bin_sequence,
                                    std::chrono::milliseconds bin_size,
                                    ctr::PcapDataBinCache* cache) {
  std::vector<ctr::TrimmedPcapDataTraceBin> bins =
      bin_sequence.AccumulateBins(bin_size, cache);
  double bins_per_second = 1000.0 / bin_size.count();

  std::vector<double> out;
  for (const auto& bin : bins) {
    double Mbits = bin.bytes * 8.0 / 1000.0 / 1000.0;
    out.emplace_back(Mbits * bins_per_second);
  }
  return out;
}

static void PlotTrace(const std::string& trace_id,
                      const ctr::BinSequence& bin_sequence,
                      ctr::PcapDataBinCache* cache) {
  std::vector<double> per_second_bits =
      BinTrace(bin_sequence, milliseconds(1000), cache);

  nc::viz::DataSeries2D to_plot;
  for (size_t i = 0; i < per_second_bits.size(); ++i) {
    to_plot.data.emplace_back(i, per_second_bits[i]);
  }
  to_plot.label = trace_id;

  nc::viz::LinePlot plot;
  plot.AddData(to_plot);
  plot.PlotToDir(nc::StrCat("binned_trace_", trace_id));

  //  std::vector<double> per_ms_bins =
  //      BinTrace(bin_sequence, milliseconds(10), cache);
  //  to_plot.data = AllanVariance(per_ms_bins);
  //
  //  nc::viz::PythonGrapher av_grapher(nc::StrCat("av_trace_", trace_id));
  //  av_grapher.PlotLine({}, {to_plot});
}

static std::vector<nc::SummaryStats> ProcessTrace(
    const ctr::PcapDataTrace* trace) {
  ctr::BinSequence bin_sequence(trace->TracesAndSlices(trace->AllSlices()));
  ctr::PcapDataBinCache cache;

  std::vector<double> per_ms = BinTrace(bin_sequence, milliseconds(1), &cache);
  std::vector<nc::SummaryStats> per_minute_stats;
  for (size_t i = 0; i < per_ms.size(); ++i) {
    size_t stats_index = i / 60000;
    per_minute_stats.resize(stats_index + 1);
    per_minute_stats[stats_index].Add(per_ms[i]);
  }

  std::string trace_id = trace->id().ToString();
  PlotTrace(trace_id, bin_sequence, &cache);
  return per_minute_stats;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  CHECK(!FLAGS_pcap_trace_store.empty());
  ctr::PcapTraceStore trace_store(FLAGS_pcap_trace_store);

  std::vector<std::string> traces_v = nc::Split(FLAGS_traces, ",");
  std::set<std::string> traces_to_plot(traces_v.begin(), traces_v.end());

  std::vector<nc::viz::DataSeries1D> deltas_to_plot;
  std::vector<nc::viz::DataSeries2D> vars_to_plot;
  for (const ctr::PcapDataTrace* trace : trace_store.AllTraces()) {
    std::string trace_id = trace->id().ToString();
    if (!traces_to_plot.empty() && !nc::ContainsKey(traces_to_plot, trace_id)) {
      continue;
    }

    LOG(INFO) << "Processing " << trace_id;
    std::vector<nc::SummaryStats> stats = ProcessTrace(trace);

    std::vector<double> deltas;
    for (size_t i = 0; i < stats.size() - 2; ++i) {
      double prev = stats[i].sum();
      double next = stats[i + 1].sum();
      deltas.emplace_back((next - prev) / prev);
    }
    nc::viz::DataSeries1D delta_data_series;
    delta_data_series.data = std::move(deltas);
    delta_data_series.label = trace_id;
    deltas_to_plot.emplace_back(delta_data_series);

    std::vector<std::pair<double, double>> vars;
    for (size_t i = 0; i < stats.size() - 2; ++i) {
      double prev_var = stats[i].var();
      double next_var = stats[i + 1].var();
      vars.emplace_back(prev_var, next_var);
    }
    nc::viz::DataSeries2D var_data_series;
    var_data_series.data = std::move(vars);
    var_data_series.label = trace_id;
    vars_to_plot.emplace_back(var_data_series);
  }

  nc::viz::CDFPlot deltas_plot;
  deltas_plot.AddData(deltas_to_plot);
  deltas_plot.PlotToDir("pcap_data_plot_deltas");

  nc::viz::LinePlot vars_plot;
  vars_plot.AddData(vars_to_plot);
  vars_plot.PlotToDir("pcap_data_plot_vars");

  return 0;
}
