#include <gflags/gflags.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/map_util.h>
#include <ncode/net/algorithm.h>
#include <ncode/net/net_common.h>
#include <ncode/strutil.h>
#include <ncode/viz/ctemplate/template.h>
#include <ncode/viz/ctemplate/template_dictionary.h>
#include <ncode/viz/ctemplate/template_enums.h>
#include <ncode/viz/grapher.h>
#include <stddef.h>
#include <stdint.h>
#include <iostream>
#include <map>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "demand_matrix_input.h"
#include "plot_algorithm_eval_tools.h"
#include "topology_input.h"

DEFINE_string(
    optimizer_labels,
    "MinMax (10 paths),LDR (no flow counts),AB4,MinMax (low delay),LDR",
    "Labels of optimizers to plot.");
DEFINE_string(output_prefix, "",
              "Prefix that will be added to all output directories.");
DECLARE_string(optimizers);

static constexpr char kTopLevelTemplate[] = "../rst/top_level_template.rst";
static constexpr char kSummaryTemplate[] = "../rst/summary_template.rst";
static constexpr char kSingleTMTemplate[] = "../rst/tm_stat_template.rst";

using namespace ctr;
using namespace std::chrono;
using namespace ctr::alg_eval;

using TopologyAndTM = std::pair<std::string, std::string>;

static std::string Indent(const std::string& input) {
  return nc::StrCat("   ", nc::StringReplace(input, "\n", "\n   ", true));
}

// Returns a list of optimizer name and optimizer label in plot.
static std::vector<std::pair<std::string, std::string>> GetOrderedOptimizers() {
  std::vector<std::string> opts = nc::Split(FLAGS_optimizers, ",");
  std::vector<std::string> labels = nc::Split(FLAGS_optimizer_labels, ",");
  CHECK(opts.size() == labels.size());

  std::vector<std::pair<std::string, std::string>> out;
  for (size_t i = 0; i < opts.size(); ++i) {
    out.emplace_back(opts[i], labels[i]);
  }
  return out;
}

class DataPlotter {
 public:
  void AddData(std::unique_ptr<TMSummary> tm_summary) {
    const nc::net::GraphStorage* graph;
    double seed;
    double load;
    double locality;
    std::tie(graph, seed, load, locality) = tm_summary->key;

    auto ll = std::make_pair(load, locality);
    data_[ll].emplace_back(tm_summary.get());
    data_storage_.emplace_back(std::move(tm_summary));
  }

  void PlotRoot(const std::string& root,
                const std::map<const nc::net::GraphStorage*,
                               nc::net::AllPairShortestPath>& sp_map) {
    using namespace std::chrono;
    ctemplate::TemplateDictionary dict("plot");

    for (const auto& load_and_locality_and_data : data_) {
      double load;
      double locality;
      std::tie(load, locality) = load_and_locality_and_data.first;

      ctemplate::TemplateDictionary* subdict =
          dict.AddSectionDictionary("subdirs");
      subdict->SetValue("name", nc::StrCat("data_load_", load, "_", locality));

      std::string subdir = nc::StrCat(root, "/data_load_", load, "_", locality);
      PlotDataByLoadAndLocality(subdir, load, locality);
    }

    for (uint32_t i = 0; i < interesting_data_.size(); ++i) {
      const TMSummary* tm = interesting_data_[i];
      const TMKey& key = tm->key;

      const nc::net::GraphStorage* graph;
      double seed;
      double load;
      double locality;
      std::tie(graph, seed, load, locality) = key;

      const nc::net::AllPairShortestPath& sp = nc::FindOrDie(sp_map, graph);
      std::string subdir = nc::StrCat(root, "/interesting_", i);
      PlotSingleTM(subdir, tm, sp, i);
    }

    std::string output;
    ctemplate::ExpandTemplate(kTopLevelTemplate, ctemplate::DO_NOT_STRIP, &dict,
                              &output);
    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root, "/index.rst"));
  }

 private:
  using LoadAndLocality = std::pair<double, double>;
  using RCSummaryVector = std::vector<std::unique_ptr<RCSummary>>;

  // Plots a single TM. The first argument is summaries, grouped by optimizer.
  // All summaries should be for the same TM.
  void PlotSingleTM(const std::string& root, const TMSummary* tm,
                    const nc::net::AllPairShortestPath& sp,
                    uint32_t interesting_index) {
    using namespace std::chrono;

    const nc::net::GraphStorage* graph;
    double seed;
    double load;
    double locality;
    std::tie(graph, seed, load, locality) = tm->key;

    const nc::lp::DemandMatrix& demand_matrix = *(tm->demand_matrix);
    SingleRCSummaryPlotPack plot_pack;
    plot_pack.PlotCumulativeDistances(demand_matrix, sp);
    plot_pack.PlotDemandSizes(demand_matrix);
    plot_pack.PlotSPUtilizations(demand_matrix);
    plot_pack.PlotTotalDelayAtLinkScale(demand_matrix);

    auto ordered_opt = GetOrderedOptimizers();
    for (const auto& opt_and_label : ordered_opt) {
      const RCSummary* rc = nc::FindOrDie(tm->rcs, opt_and_label.first).get();
      plot_pack.PlotDelayVsUtilization(opt_and_label.first, *(rc));
      plot_pack.PlotAbsoluteStretch(opt_and_label.first, *(rc));
    }

    ctemplate::TemplateDictionary dict("plot");
    PlotAndAddToTemplate(root, "demand_sizes", plot_pack.demand_sizes_plot(),
                         &dict);
    PlotAndAddToTemplate(root, "sp_utilization",
                         plot_pack.sp_utilizations_plot(), &dict);
    PlotAndAddToTemplate(root, "total_delay_at_link_scale",
                         plot_pack.total_delay_at_link_scale_plot(), &dict);
    PlotAndAddToTemplate(root, "cumulative_demands",
                         plot_pack.cumulative_demands_plot(), &dict);
    PlotAndAddToTemplate(root, "cumulative_demands_hop",
                         plot_pack.cumulative_demands_hop_plot(), &dict);
    PlotAndAddToTemplate(root, "link_delay_vs_utilization",
                         plot_pack.link_delay_vs_link_utilization_plot(),
                         &dict);
    PlotAndAddToTemplate(root, "absolute_stretch",
                         plot_pack.absolute_path_stretch_plot(), &dict);

    dict.SetValue(
        "demand_sizes_interesting_table",
        plot_pack.demand_sizes_interesting().ToRSTTable(&interesting_data_));

    std::string name_stripped =
        nc::File::ExtractFileName(*(tm->topology_file_name));
    name_stripped = nc::StringReplace(name_stripped, ".graph", "", true);

    dict.SetValue("topology_name", name_stripped);
    dict.SetValue("tm_seed", nc::StrCat(seed));
    dict.SetValue("load", nc::StrCat(load));
    dict.SetValue("locality", nc::StrCat(locality));
    dict.SetValue("interesting_index", nc::StrCat(interesting_index));
    dict.SetValue("tm_summary", Indent(demand_matrix.ToString()));

    std::string output;
    ctemplate::ExpandTemplate(kSingleTMTemplate, ctemplate::DO_NOT_STRIP, &dict,
                              &output);
    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root, "/index.rst"));
  }

  void PlotAndAddToTemplate(const std::string& root, const std::string& base,
                            const nc::viz::Plot& plot,
                            ctemplate::TemplateDictionary* dict) {
    std::string plot_location = nc::StrCat(root, "/", base, ".svg");
    std::string tgz_location = nc::StrCat(root, "/", base, ".tgz");
    plot.PlotToSVGFile(plot_location);
    plot.PlotToArchiveFile(tgz_location);

    std::string rel_location = nc::StrCat(base, ".svg");
    std::string rel_tgz_location = nc::StrCat(base, ".tgz");
    dict->SetValue(nc::StrCat(base, "_location"), rel_location);
    dict->SetValue(nc::StrCat(base, "_location_tgz"), rel_tgz_location);
  }

  std::map<std::string, std::vector<const RCSummary*>> DataByOptimizer(
      double load, double locality) {
    const std::vector<const TMSummary*>& all_data =
        nc::FindOrDieNoPrint(data_, std::make_pair(load, locality));

    std::map<std::string, std::vector<const RCSummary*>> out;
    for (const TMSummary* tm : all_data) {
      for (const auto& opt_and_rc : tm->rcs) {
        const std::string& opt = opt_and_rc.first;
        const RCSummary* rc = opt_and_rc.second.get();
        out[opt].emplace_back(rc);
      }
    }
    return out;
  }

  void PlotDataByLoadAndLocality(const std::string& root, double load,
                                 double locality) {
    MultiRCSummaryPlotPack plot_pack;
    ctemplate::TemplateDictionary dict("plot");
    std::map<std::string, std::vector<const RCSummary*>> opt_map =
        DataByOptimizer(load, locality);

    auto ordered_opt = GetOrderedOptimizers();
    for (const auto& opt_and_label : ordered_opt) {
      const std::vector<const RCSummary*>& rcs =
          nc::FindOrDie(opt_map, opt_and_label.first);
      const std::string& label = opt_and_label.second;
      plot_pack.PlotLinkUtilizations(label, rcs);
      plot_pack.PlotPathCounts(label, rcs);
      plot_pack.PlotPathRatios(label, rcs);
      plot_pack.PlotStretchDistribution(label, rcs);
    }

    const std::vector<const TMSummary*>& tms =
        nc::FindOrDieNoPrint(data_, std::make_pair(load, locality));
    plot_pack.PlotLinkScalesAtDelay(tms);

    PlotAndAddToTemplate(root, "path_ratios", plot_pack.ratios_plot(), &dict);
    PlotAndAddToTemplate(root, "path_stretch_rel",
                         plot_pack.path_stretch_rel_plot(), &dict);
    PlotAndAddToTemplate(root, "max_path_stretch_rel",
                         plot_pack.path_stretch_max_rel_plot(), &dict);
    PlotAndAddToTemplate(root, "path_count",
                         plot_pack.aggregate_path_count_plot(), &dict);
    PlotAndAddToTemplate(root, "link_utilization",
                         plot_pack.link_utilization_plot(), &dict);
    PlotAndAddToTemplate(root, "link_scales", plot_pack.link_scales_plot(),
                         &dict);

    dict.SetValue(
        "path_ratios_interesting_table",
        plot_pack.ratios_interesting().ToRSTTable(&interesting_data_));
    dict.SetValue("max_path_stretch_rel_interesting_table",
                  plot_pack.path_stretch_max_rel_interesting().ToRSTTable(
                      &interesting_data_));
    dict.SetValue(
        "link_scales_interesting_table",
        plot_pack.link_scales_interesting().ToRSTTable(&interesting_data_));

    dict.SetValue("tm_load", nc::StrCat(load));
    dict.SetValue("tm_locality", nc::StrCat(locality));
    std::string output;
    ctemplate::ExpandTemplate(kSummaryTemplate, ctemplate::DO_NOT_STRIP, &dict,
                              &output);
    nc::File::WriteStringToFileOrDie(output, nc::StrCat(root, "/index.rst"));
  }

  // All data, grouped by load and locality.
  std::map<LoadAndLocality, std::vector<const TMSummary*>> data_;

  std::vector<const TMSummary*> interesting_data_;

  // All data is stored here.
  std::vector<std::unique_ptr<TMSummary>> data_storage_;
};

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::vector<std::unique_ptr<TMSummary>> summaries;

  std::tie(topologies, matrices, summaries) = GetTMSummaries();

  // Need to know the shortest paths for each graph.
  std::map<const nc::net::GraphStorage*, nc::net::AllPairShortestPath> sp_map;

  std::map<std::string, const ctr::TopologyAndFilename*> topologies_by_name;
  for (const auto& topology : topologies) {
    topologies_by_name[topology.file] = &topology;

    const nc::net::GraphStorage* graph = topology.graph.get();
    sp_map.emplace(
        std::piecewise_construct, std::forward_as_tuple(graph),
        std::forward_as_tuple(nc::net::ExclusionSet(), graph->AdjacencyList(),
                              nullptr, nullptr));
  }

  // A tree of plots.
  DataPlotter data_plotter;
  for (auto& summary : summaries) {
    data_plotter.AddData(std::move(summary));
  }

  data_plotter.PlotRoot("plot_tree", sp_map);
}
