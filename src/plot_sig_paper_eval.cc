#include <gflags/gflags.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/map_util.h>
#include <ncode/net/net_common.h>
#include <ncode/strutil.h>
#include <ncode/viz/grapher.h>
#include <stddef.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "demand_matrix_input.h"
#include "opt/opt.h"
#include "plot_algorithm_eval_tools.h"
#include "topology_input.h"

DEFINE_string(output_prefix, "sig",
              "Prefix that will be added to all output directories.");
DEFINE_double(
    routability_threshold, 0.5,
    "Threshold for routability, will do separate plots for topologies "
    "below and above that");
DEFINE_double(sp_fraction, 1.4, "How far from the SP a path can be");
DEFINE_double(link_fraction_limit, 0.7,
              "At least this much of the SP's links can be routed around");
DEFINE_bool(dump_median_max_stretch, false,
            "If true will not plot anything, but will only dump the median max "
            "path stretch");

using namespace ctr;
using namespace std::chrono;
using namespace ctr::alg_eval;

// Groups a series of TMSummary by optimizer. Will only consider graphs in the
// given set.
std::map<std::string, std::vector<const RCSummary*>> GroupByOpt(
    std::set<const nc::net::GraphStorage*>& graphs,
    const std::vector<std::unique_ptr<TMSummary>>& tm_summaries) {
  std::map<std::string, std::vector<const RCSummary*>> out;
  for (const auto& tm_summary : tm_summaries) {
    const nc::net::GraphStorage* graph = tm_summary->demand_matrix->graph();
    if (!nc::ContainsKey(graphs, graph)) {
      continue;
    }

    for (const auto& opt_and_summary : tm_summary->rcs) {
      out[opt_and_summary.first].emplace_back(opt_and_summary.second.get());
    }
  }

  return out;
}

void PlotStuff(
    const std::string& prefix,
    const std::map<std::string, std::vector<const RCSummary*>>& summaries) {
  std::map<size_t, nc::viz::CDFPlot> stretch_plots;
  nc::viz::CDFPlot total_delay_delta_plot(
      {"Distribution of total delay", "change over SP delay"});
  nc::viz::CDFPlot avg_path_count_plot(
      {"Average paths per aggregate", "average path count"});
  nc::viz::CDFPlot max_path_count_plot(
      {"Max paths per aggregate", "max path count"});
  nc::viz::CDFPlot single_path_fraction_plot(
      {"Fraction of aggregates with a single path", "fraction"});
  nc::viz::CDFPlot sp_fraction_plot(
      {"Fraction of aggregates with stretch of 1", "fraction"});

  std::vector<uint32_t> percentiles = {90, 95, 100};
  for (uint32_t i : percentiles) {
    stretch_plots[i] = nc::viz::CDFPlot(
        {nc::Substitute("Distribution of the $0th percentile of flow stretch",
                        i),
         "times greater than shortest path"});
  }

  for (const auto& opt_and_summaries : summaries) {
    const std::string& opt = opt_and_summaries.first;
    const std::vector<const RCSummary*>& rcs = opt_and_summaries.second;

    std::vector<std::pair<double, const RCSummary*>> delay_changes;
    std::vector<double> avg_path_counts;
    std::vector<double> max_path_counts;
    std::vector<double> fractions_on_sp;
    std::vector<double> fractions_on_single_path;
    std::map<size_t, std::vector<double>> values;
    for (const RCSummary* rc : rcs) {
      std::vector<double> stretches = rc->rel_stretches;
      std::vector<bool> overloaded = rc->overloaded;
      std::sort(stretches.begin(), stretches.end());

      double all_flows =
          std::accumulate(rc->flow_counts.begin(), rc->flow_counts.end(), 0.0);
      double flows_on_sp = 0;
      for (size_t i = 0; i < rc->abs_stretches.size(); ++i) {
        if (rc->abs_stretches[i] < 0.0001) {
          flows_on_sp += rc->flow_counts[i];
        }
      }
      fractions_on_sp.emplace_back(flows_on_sp / all_flows);

      double aggregates_on_single_path = 0;
      for (size_t i = 0; i < rc->path_counts.size(); ++i) {
        if (rc->path_counts[i] == 1) {
          ++aggregates_on_single_path;
        }
      }
      fractions_on_single_path.emplace_back(aggregates_on_single_path /
                                            rc->path_counts.size());

      double total =
          std::accumulate(rc->path_counts.begin(), rc->path_counts.end(), 0.0);
      double max_path_count =
          *(std::max_element(rc->path_counts.begin(), rc->path_counts.end()));

      avg_path_counts.emplace_back(total / rc->path_counts.size());
      max_path_counts.emplace_back(max_path_count);

      bool any_overloaded = false;
      CHECK(stretches.size() == overloaded.size());
      for (size_t i = 0; i < overloaded.size(); ++i) {
        if (overloaded[i]) {
          any_overloaded = true;
          break;
        }
      }

      if (any_overloaded) {
        continue;
      }

      std::vector<double> p = nc::Percentiles(&stretches);
      for (size_t i : percentiles) {
        values[i].emplace_back(p[i]);
      }

      double fraction_change_in_delay =
          (rc->total_delay - rc->total_sp_delay) / rc->total_delay;
      delay_changes.emplace_back(fraction_change_in_delay, rc);
    }

    std::sort(delay_changes.begin(), delay_changes.end());
    const RCSummary* median_delay_change =
        delay_changes[delay_changes.size() / 2].second;

    std::vector<double> delay_change_values;
    for (const auto& delay_change : delay_changes) {
      delay_change_values.emplace_back(delay_change.first);
    }

    LOG(INFO) << "Median for " << opt << " is "
              << (*median_delay_change->parent->topology_file_name) << " tm "
              << (*median_delay_change->parent->demand_matrix_file_name);

    total_delay_delta_plot.AddData(opt, delay_change_values);
    avg_path_count_plot.AddData(opt, avg_path_counts);
    max_path_count_plot.AddData(opt, max_path_counts);
    sp_fraction_plot.AddData(opt, fractions_on_sp);
    single_path_fraction_plot.AddData(opt, fractions_on_single_path);
    for (const auto& percentile_and_values : values) {
      size_t percentile = percentile_and_values.first;
      const std::vector<double>& values = percentile_and_values.second;
      stretch_plots[percentile].AddData(opt, values);
    }
  }

  for (const auto& percentile_and_plot : stretch_plots) {
    uint32_t percentile = percentile_and_plot.first;
    const nc::viz::CDFPlot& plot = percentile_and_plot.second;
    plot.PlotToDir(nc::StrCat(prefix, "p", percentile, "_stretch"));
  }

  total_delay_delta_plot.PlotToDir(nc::StrCat(prefix, "total_delay_delta"));
  avg_path_count_plot.PlotToDir(nc::StrCat(prefix, "avg_path_count"));
  max_path_count_plot.PlotToDir(nc::StrCat(prefix, "max_path_count"));
  sp_fraction_plot.PlotToDir(nc::StrCat(prefix, "sp_fraction"));
  single_path_fraction_plot.PlotToDir(
      nc::StrCat(prefix, "single_path_fraction"));
}

static void PlotRuntime(
    const std::string& prefix,
    const std::map<std::string, std::vector<const RCSummary*>>& summaries) {
  nc::viz::CDFPlot runtime_plot = nc::viz::CDFPlot({"Runtime", "time (ms)"});

  for (const auto& opt_and_summaries : summaries) {
    const std::string& opt = opt_and_summaries.first;
    const std::vector<const RCSummary*>& rcs = opt_and_summaries.second;

    std::vector<double> values;
    for (const RCSummary* rc : rcs) {
      values.emplace_back(rc->runtime_ms);
    }
    runtime_plot.AddData(opt, values);
  }

  runtime_plot.PlotToDir(nc::StrCat(prefix, "runtime"));
}

static double MedianOfPercentileStretch(
    const std::vector<const RCSummary*>& rcs, size_t p) {
  CHECK(p < 101);
  std::vector<double> values;
  for (const RCSummary* rc : rcs) {
    std::vector<double> stretches = rc->rel_stretches;
    std::vector<bool> overloaded = rc->overloaded;
    std::sort(stretches.begin(), stretches.end());

    bool any_overloaded = false;
    CHECK(stretches.size() == overloaded.size());
    for (size_t i = 0; i < overloaded.size(); ++i) {
      if (overloaded[i]) {
        any_overloaded = true;
        break;
      }
    }

    if (any_overloaded) {
      values.emplace_back(std::numeric_limits<double>::max());
      continue;
    }

    std::vector<double> stretches_p = nc::Percentiles(&stretches);
    values.emplace_back(stretches_p[p]);
  }

  std::sort(values.begin(), values.end());
  CHECK(!values.empty());
  return values[values.size() / 2];
}

static void DumpMeidanOfPercentileStretch(
    const std::map<std::string, std::vector<const RCSummary*>>& summaries,
    size_t p) {
  std::vector<std::string> medians;
  for (const auto& opt_and_summaries : summaries) {
    const std::string& opt = opt_and_summaries.first;
    const std::vector<const RCSummary*>& summaries_for_opt =
        opt_and_summaries.second;
    double median = MedianOfPercentileStretch(summaries_for_opt, p);
    medians.emplace_back(nc::StrCat("('", opt, "',", median, ")"));
  }

  std::string out = nc::StrCat("[", nc::Join(medians, ","), "]");
  nc::File::WriteStringToFileOrDie(out, "stretch");
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::vector<std::unique_ptr<TMSummary>> summaries;
  std::tie(topologies, matrices, summaries) = GetTMSummaries();

  // Need to know the routability paths for each graph.
  std::set<const nc::net::GraphStorage*> high_routability_graphs;
  std::set<const nc::net::GraphStorage*> low_routability_graphs;
  for (const auto& topology : topologies) {
    const nc::net::GraphStorage* graph = topology.graph.get();
    double routability = GetFractionOfPairsAboveLinkFraction(
        *graph, FLAGS_sp_fraction, FLAGS_link_fraction_limit);

    if (routability > FLAGS_routability_threshold) {
      high_routability_graphs.insert(graph);
    } else {
      low_routability_graphs.insert(graph);
    }
  }

  auto hr_by_opt = GroupByOpt(high_routability_graphs, summaries);
  auto lr_by_opt = GroupByOpt(low_routability_graphs, summaries);

  if (FLAGS_dump_median_max_stretch) {
    DumpMeidanOfPercentileStretch(hr_by_opt, 100);
    return 0;
  }

  std::string hr_prefix = nc::StrCat(FLAGS_output_prefix, "_hr_");
  std::string lr_prefix = nc::StrCat(FLAGS_output_prefix, "_lr_");
  PlotStuff(hr_prefix, hr_by_opt);
  PlotStuff(lr_prefix, lr_by_opt);
  PlotRuntime(hr_prefix, hr_by_opt);
  PlotRuntime(lr_prefix, lr_by_opt);
}
