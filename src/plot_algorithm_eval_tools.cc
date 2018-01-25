#include "plot_algorithm_eval_tools.h"

#include <gflags/gflags.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/map_util.h>
#include <ncode/net/algorithm.h>
#include <ncode/perfect_hash.h>
#include <ncode/stats.h>
#include <ncode/strutil.h>
#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <type_traits>
#include <utility>

#include "common.h"
#include "opt/path_provider.h"

DEFINE_string(optimizers, "MinMaxK10,CTRNFC,B4,MinMaxLD,CTR",
              "Optimizers to parse.");
DEFINE_bool(skip_trivial, false, "Whether or not to skip trivial matrices.");
DEFINE_uint64(delay_penalty_ms, 10,
              "Delay penalty added to all paths that cross congested links.");

static constexpr size_t kDiscreteMultiplier = 10000;
static constexpr size_t kPercentilesCount = 10000;

using namespace ctr;
using namespace std::chrono;

namespace ctr {
namespace alg_eval {

static double ToMs(nc::net::Delay delay) {
  return duration_cast<microseconds>(delay).count() / 1000.0;
}

static std::string GetFilename(const std::string& tm_file,
                               const std::string opt_string) {
  std::string tm_base = nc::StringReplace(tm_file, ".demands", "", true);
  return nc::StrCat(tm_base, "_", opt_string, ".rc");
}

static std::unique_ptr<TMSummary> ParseTM(
    const nc::lp::DemandMatrix* demand_matrix, const std::string* top_file,
    const std::string* tm_file) {
  auto out = nc::make_unique<TMSummary>();
  double locality;
  std::string locality_str =
      nc::FindOrDie(demand_matrix->properties(), "locality");
  CHECK(nc::safe_strtod(locality_str, &locality));

  double load;
  std::string load_str = nc::FindOrDie(demand_matrix->properties(), "load");
  CHECK(nc::safe_strtod(load_str, &load));

  double seed;
  std::string seed_str = nc::FindOrDie(demand_matrix->properties(), "seed");
  CHECK(nc::safe_strtod(seed_str, &seed));

  out->demand_matrix = demand_matrix;
  out->key = std::make_tuple(demand_matrix->graph(), seed, load, locality);
  out->topology_file_name = top_file;
  out->demand_matrix_file_name = tm_file;
  return out;
}

static std::unique_ptr<RCSummary> ParseRC(
    const RoutingConfiguration& rc, const nc::net::AllPairShortestPath& sp,
    const std::string* opt) {
  auto out = nc::make_unique<RCSummary>();

  nc::net::GraphLinkMap<double> utilizations = rc.LinkUtilizations();
  out->utilizations = utilizations;

  nc::net::GraphLinkSet overloaded_links;
  for (const auto& link_and_utilization : utilizations) {
    nc::net::GraphLinkIndex link = link_and_utilization.first;
    double utilization = *(link_and_utilization.second);
    if (utilization > 1.001) {
      overloaded_links.insert(link);
    }
  }

  for (const auto& aggregate_and_routes : rc.routes()) {
    const AggregateId& id = aggregate_and_routes.first;
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(rc.demands(), id);
    double aggregate_flow_count = demand_and_flow_count.second;

    nc::net::Delay sp_delay = sp.GetDistance(id.src(), id.dst());
    double sp_delay_ms = ToMs(sp_delay);
    sp_delay_ms = std::max(1.0, sp_delay_ms);

    for (const RouteAndFraction& route_and_fraction :
         aggregate_and_routes.second) {
      const nc::net::Walk* path = route_and_fraction.first;
      double fraction = route_and_fraction.second;

      double path_delay_ms = ToMs(path->delay());
      path_delay_ms = std::max(1.0, path_delay_ms);
      bool overloaded = path->ContainsAny(overloaded_links);
      if (overloaded) {
        path_delay_ms += FLAGS_delay_penalty_ms;
      }

      double abs_stretch = path_delay_ms - sp_delay_ms;
      double rel_stretch = path_delay_ms / sp_delay_ms;
      CHECK(abs_stretch >= 0);
      CHECK(rel_stretch >= 0);
      CHECK(rel_stretch < 1000);

      size_t flow_count = aggregate_flow_count * fraction;
      flow_count = std::max(1ul, flow_count);
      out->abs_stretches.emplace_back(abs_stretch);
      out->rel_stretches.emplace_back(rel_stretch);
      out->flow_counts.emplace_back(flow_count);
      out->total_sp_delay += sp_delay_ms * flow_count;
      out->total_delay += (sp_delay_ms + abs_stretch) * flow_count;
    }
    out->path_counts.emplace_back(aggregate_and_routes.second.size());
  }

  out->opt = opt;
  return out;
}

std::tuple<std::vector<TopologyAndFilename>,
           std::vector<DemandMatrixAndFilename>,
           std::vector<std::unique_ptr<TMSummary>>>
GetTMSummaries() {
  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::tie(topologies, matrices) =
      ctr::GetDemandMatrixInputs(FLAGS_skip_trivial);
  std::vector<std::string> optimizers = nc::Split(FLAGS_optimizers, ",");

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

  std::vector<std::unique_ptr<TMSummary>> summaries;
  for (const auto& matrix : matrices) {
    const ctr::TopologyAndFilename* topology_and_filename =
        nc::FindOrDieNoPrint(topologies_by_name, matrix.topology_file);

    ctr::PathProvider path_provider(topology_and_filename->graph.get());
    const nc::lp::DemandMatrix* demand_matrix = matrix.demand_matrix.get();

    std::vector<std::unique_ptr<RoutingConfiguration>> rcs;
    for (const std::string& opt : optimizers) {
      std::string rc_filename = GetFilename(matrix.file, opt);
      if (!nc::File::Exists(rc_filename)) {
        LOG(INFO) << "Missing " << rc_filename;
        break;
      }

      std::string rc_serialized = nc::File::ReadFileToStringOrDie(rc_filename);
      std::unique_ptr<TrafficMatrix> tm =
          TrafficMatrix::DistributeFromDemandMatrix(*demand_matrix);
      LOG(INFO) << "Will parse " << matrix.file << " at "
                << matrix.topology_file << " opt " << opt;
      auto rc = RoutingConfiguration::LoadFromSerializedText(
          *tm, topology_and_filename->node_order, rc_serialized,
          &path_provider);
      rcs.emplace_back(std::move(rc));
    }

    if (rcs.size() != optimizers.size()) {
      LOG(INFO) << "Will skip " << matrix.file << " at "
                << matrix.topology_file;
      continue;
    }

    auto tm_summary =
        ParseTM(demand_matrix, &(matrix.file), &(matrix.topology_file));
    for (size_t i = 0; i < optimizers.size(); ++i) {
      const std::string& opt = optimizers[i];
      const RoutingConfiguration& rc = *(rcs[i]);

      const nc::net::AllPairShortestPath& sp =
          nc::FindOrDieNoPrint(sp_map, rc.graph());
      auto summary = ParseRC(rc, sp, &opt);
      summary->parent = tm_summary.get();
      tm_summary->rcs[opt] = std::move(summary);
    }
    summaries.emplace_back(std::move(tm_summary));
  }

  return std::make_tuple(std::move(topologies), std::move(matrices),
                         std::move(summaries));
}

static std::vector<std::string> Escape(const std::vector<std::string>& row) {
  std::vector<std::string> new_row;
  for (const auto& element : row) {
    new_row.emplace_back(nc::StrCat("\"", element, "\""));
  }
  return new_row;
}

void InterestingTable::Add(const std::vector<std::string>& row,
                           const TMSummary* tm) {
  rows_.emplace_back(Escape(row));
  tms_.emplace_back(tm);
}

void InterestingTable::Add(const std::vector<std::string>& row) {
  rows_.emplace_back(Escape(row));
}

std::string InterestingTable::ToRSTTable(
    std::vector<const TMSummary*>* interestring_tm_list) const {
  bool has_rcs = !tms_.empty();
  if (has_rcs) {
    CHECK(rows_.size() == tms_.size());
  }

  std::string out = nc::StrCat(".. csv-table:: ", title_, "\n");
  nc::StrAppend(&out, "   :header: ", nc::Join(header_, ", "));
  if (has_rcs) {
    nc::StrAppend(&out, ", link");
  }
  nc::StrAppend(&out, "\n");

  std::vector<std::string> widths;
  size_t count = has_rcs ? header_.size() + 1 : header_.size();
  for (size_t i = 0; i < count; ++i) {
    widths.emplace_back("1");
  }

  nc::StrAppend(&out, "   :widths: ", nc::Join(widths, ","), "\n\n");
  for (size_t row_i = 0; row_i < rows_.size(); ++row_i) {
    const std::vector<std::string>& row = rows_[row_i];
    nc::StrAppend(&out, "   ", nc::Join(row, ","));

    if (has_rcs) {
      uint32_t current_size = interestring_tm_list->size();
      nc::StrAppend(&out, ", see :ref:`interesting_", current_size, "`");
      interestring_tm_list->emplace_back(tms_[row_i]);
    }
    nc::StrAppend(&out, "\n");
  }

  return out;
}

// Parses a comma-separated string of doubles.
static std::vector<double> GetCommaSeparated(const std::string& value) {
  std::vector<std::string> split = nc::Split(value, ",");
  std::vector<double> out;
  for (const auto& piece : split) {
    double v;
    CHECK(nc::safe_strtod(piece, &v));
    out.emplace_back(v);
  }

  return out;
}

// Returns the link scale at which total delay will be 'fraction' of the one at
// no headroom.
static double LinkScaleAtDelayFraction(
    double fraction, const nc::lp::DemandMatrix& demand_matrix) {
  const std::map<std::string, std::string>& properties =
      demand_matrix.properties();
  std::string step_str = nc::FindOrDie(properties, "headroom_vs_delay_step");
  std::string delays_at_scale_str =
      nc::FindOrDie(properties, "headroom_vs_delay");

  double abs_stride;
  CHECK(nc::safe_strtod(step_str, &abs_stride));
  const std::vector<double>& delay_at_scale =
      GetCommaSeparated(delays_at_scale_str);

  std::vector<double> xs;
  std::vector<double> ys;
  for (size_t i = 0; i < delay_at_scale.size(); ++i) {
    double multiplier = 1.0 - i * abs_stride;
    xs.emplace_back(multiplier);
    ys.emplace_back(delay_at_scale[i]);
  }

  nc::Empirical2DFunction f(xs, ys, nc::Empirical2DFunction::LINEAR);
  double min_delay = delay_at_scale.front();
  double max_delay = delay_at_scale.back();
  if (max_delay / min_delay <= fraction) {
    return 1.0 - delay_at_scale.size() * abs_stride;
  }

  for (double link_scale = 1.0; link_scale > 0.0; link_scale -= 0.01) {
    double delay = f.Eval(link_scale);
    double delay_fraction = delay / min_delay;
    if (delay_fraction >= fraction) {
      return link_scale;
    }
  }
  LOG(FATAL) << "Should not happen";
  return 0;
}

MultiRCSummaryPlotPack::MultiRCSummaryPlotPack()
    : path_stretch_rel_plot_(
          {"Relative path stretch", "fraction longer than SP", "CDF"}),
      path_stretch_max_rel_plot_(
          {"Maximum relative path stretch", "fraction longer than SP"}),
      path_stretch_max_rel_interesting_({"optimizer", "type", "value"}),
      link_utilization_plot_({"Link utilizations", "utilization"}),
      link_utilization_interesting_({"optimizer", "type", "value"}),
      ratios_plot_({"Total delay", "fraction total delay is longer than SP"}),
      ratios_interesting_({"optimizer", "type", "value"}),
      aggregate_path_count_plot_(
          {"Number of paths per aggregate", "path count"}),
      aggregate_path_count_interesting_({"optimizer", "type", "value"}),
      link_scales_plot_(
          {"Link scales at X% increase in total delay", "link scale"}),
      link_scales_interesting_({"percent", "type", "value"}) {}

void MultiRCSummaryPlotPack::PlotStretchDistribution(
    const std::string& optimizer, const std::vector<const RCSummary*>& rcs) {
  nc::DiscreteDistribution<uint64_t> dist;
  std::vector<std::pair<double, const RCSummary*>> maxs_decorated;

  for (const RCSummary* rc : rcs) {
    double max = 0;
    for (size_t i = 0; i < rc->rel_stretches.size(); ++i) {
      double rel_stretch = rc->rel_stretches[i];
      double flow_count = rc->flow_counts[i];
      uint64_t value_discrete =
          static_cast<uint64_t>(kDiscreteMultiplier * rel_stretch);
      dist.Add(value_discrete, flow_count);
      max = std::max(rel_stretch, max);
    }
    maxs_decorated.emplace_back(max, rc);
  }

  std::sort(maxs_decorated.begin(), maxs_decorated.end());
  const auto* max = &maxs_decorated.back();
  const auto* med = &maxs_decorated[maxs_decorated.size() / 2];
  path_stretch_max_rel_interesting_.Add(
      {optimizer, "max", nc::ToStringMaxDecimals(max->first, 2)},
      max->second->parent);
  path_stretch_max_rel_interesting_.Add(
      {optimizer, "med", nc::ToStringMaxDecimals(med->first, 2)},
      med->second->parent);

  std::vector<uint64_t> percentiles = dist.Percentiles(kPercentilesCount);
  CHECK(percentiles.size() == kPercentilesCount + 1);

  std::vector<std::pair<double, double>> path_stretch_rel_data;
  for (size_t i = 0; i < percentiles.size(); ++i) {
    double discrete_value = percentiles[i];
    double p = discrete_value / kDiscreteMultiplier;
    path_stretch_rel_data.emplace_back(
        p, static_cast<double>(i) / kPercentilesCount);
  }
  path_stretch_rel_plot_.AddData(optimizer, path_stretch_rel_data);

  std::vector<double> maxs;
  for (const auto& value_and_rc : maxs_decorated) {
    maxs.emplace_back(value_and_rc.first);
  }
  path_stretch_max_rel_plot_.AddData(optimizer, maxs);
}

void MultiRCSummaryPlotPack::PlotPathRatios(
    const std::string& optimizer, const std::vector<const RCSummary*>& rcs) {
  std::vector<std::pair<double, const RCSummary*>> decorated;
  for (const RCSummary* rc : rcs) {
    double change = (rc->total_delay - rc->total_sp_delay) / rc->total_sp_delay;
    decorated.emplace_back(change, rc);
  }

  std::sort(decorated.begin(), decorated.end());
  const RCSummary* max = decorated.back().second;
  const RCSummary* med = decorated[decorated.size() / 2].second;

  std::vector<double> ratios;
  for (const auto& change_and_rc : decorated) {
    ratios.emplace_back(change_and_rc.first);
  }

  ratios_interesting_.Add(
      {optimizer, "max", nc::ToStringMaxDecimals(ratios.back(), 2)},
      max->parent);
  ratios_interesting_.Add(
      {optimizer, "med", nc::ToStringMaxDecimals(ratios[ratios.size() / 2], 2)},
      med->parent);
  ratios_plot_.AddData(optimizer, ratios);
}

void MultiRCSummaryPlotPack::PlotPathCounts(
    const std::string& optimizer, const std::vector<const RCSummary*>& rcs) {
  std::vector<std::pair<double, const RCSummary*>> decorated;
  for (const RCSummary* rc : rcs) {
    for (double path_count : rc->path_counts) {
      decorated.emplace_back(path_count, rc);
    }
  }

  std::sort(decorated.begin(), decorated.end());
  const RCSummary* max = decorated.back().second;
  const RCSummary* med = decorated[decorated.size() / 2].second;

  std::vector<double> path_counts;
  for (const auto& value_and_rc : decorated) {
    path_counts.emplace_back(value_and_rc.first);
  }

  aggregate_path_count_interesting_.Add(
      {optimizer, "max", nc::ToStringMaxDecimals(path_counts.back(), 2)},
      max->parent);
  aggregate_path_count_interesting_.Add(
      {optimizer, "med",
       nc::ToStringMaxDecimals(path_counts[path_counts.size() / 2], 2)},
      med->parent);
  aggregate_path_count_plot_.AddData(optimizer, path_counts);
}

void MultiRCSummaryPlotPack::PlotLinkUtilizations(
    const std::string& optimizer, const std::vector<const RCSummary*>& rcs) {
  std::vector<std::pair<double, const RCSummary*>> decorated;
  for (const RCSummary* rc : rcs) {
    for (const auto& link_and_utilization : rc->utilizations) {
      double link_utilization = *(link_and_utilization.second);
      decorated.emplace_back(link_utilization, rc);
    }
  }

  std::sort(decorated.begin(), decorated.end());
  const RCSummary* max = decorated.back().second;
  const RCSummary* med = decorated[decorated.size() / 2].second;

  std::vector<double> link_utilizations;
  for (const auto& value_and_rc : decorated) {
    link_utilizations.emplace_back(value_and_rc.first);
  }

  link_utilization_interesting_.Add(
      {optimizer, "max", nc::ToStringMaxDecimals(link_utilizations.back(), 2)},
      max->parent);
  link_utilization_interesting_.Add(
      {optimizer, "med",
       nc::ToStringMaxDecimals(link_utilizations[link_utilizations.size() / 2],
                               2)},
      med->parent);
  link_utilization_plot_.AddData(optimizer, link_utilizations);
}

void MultiRCSummaryPlotPack::PlotLinkScalesAtDelay(
    const std::vector<const TMSummary*>& all_demands) {
  std::vector<double> p_values = {1.01, 1.02, 1.05, 1.10};
  for (double p : p_values) {
    std::vector<std::pair<double, const TMSummary*>> decorated;

    for (const TMSummary* tm : all_demands) {
      double scale = LinkScaleAtDelayFraction(p, *(tm->demand_matrix));
      decorated.emplace_back(scale, tm);
    }

    std::sort(decorated.begin(), decorated.end());
    const TMSummary* max = decorated.back().second;
    const TMSummary* med = decorated[decorated.size() / 2].second;

    std::vector<double> link_scales;
    for (const auto& value_and_tm : decorated) {
      link_scales.emplace_back(value_and_tm.first);
    }

    uint32_t percents = (p - 1) * 100;
    std::string p_string = nc::StrCat(percents, "%");

    link_scales_interesting_.Add(
        {p_string, "max", nc::ToStringMaxDecimals(link_scales.back(), 2)}, max);
    link_scales_interesting_.Add(
        {p_string, "med",
         nc::ToStringMaxDecimals(link_scales[link_scales.size() / 2], 2)},
        med);
    link_scales_plot_.AddData(p_string, link_scales);
  }
}

SingleRCSummaryPlotPack::SingleRCSummaryPlotPack()
    : demand_sizes_plot_(
          {"Demand sizes (all possible demands considered)", "size (Mbps)"}),
      demand_sizes_interesting_({"type", "value", "demand"}),
      sp_utilizations_plot_(
          {"Shortest path link utilization", "link utilization"}),
      sp_utilizations_interesting_({"type", "utilization", "link"}),
      cumulative_demands_plot_({"Cumulative distance vs aggregate size",
                                "SP distance (ms)", "Cumulative size (Mbps)"}),
      cumulative_demands_hop_plot_({"Cumulative distance vs aggregate size",
                                    "SP distance (hops)",
                                    "Cumulative size (Mbps)"}),
      total_delay_at_link_scale_plot_(
          {"Total delay vs link scale", "Link scale", "Total delay"}),
      link_delay_vs_link_utilization_plot_({"Link delay vs link uilization",
                                            "Links ranked by propagation delay",
                                            "Link utilization"}),
      absolute_path_stretch_plot_(
          {"Absolute path stretch", "path stretch (ms)", "CDF"}) {
  link_delay_vs_link_utilization_plot_.TurnIntoScatterPlot();
}

void SingleRCSummaryPlotPack::PlotCumulativeDistances(
    const nc::lp::DemandMatrix& demand_matrix,
    const nc::net::AllPairShortestPath& sp) {
  std::vector<std::pair<double, double>> distances_and_demands;

  // Since hop counts are quantized, will store them in a map.
  std::map<uint32_t, double> hop_to_demand;

  for (const auto& element : demand_matrix.elements()) {
    nc::net::GraphNodeIndex src = element.src;
    nc::net::GraphNodeIndex dst = element.dst;
    double demand_Mbps = element.demand.Mbps();
    if (demand_Mbps == 0) {
      continue;
    }

    auto path = sp.GetPath(src, dst);
    double distance_ms = duration_cast<milliseconds>(path->delay()).count();
    distances_and_demands.emplace_back(distance_ms, demand_Mbps);

    uint32_t hop_count = path->links().size();
    hop_to_demand[hop_count] += demand_Mbps;
  }
  std::sort(distances_and_demands.begin(), distances_and_demands.end());

  double total = 0;
  for (auto& distance_and_demand : distances_and_demands) {
    double& demand = distance_and_demand.second;
    total += demand;
    demand = total;
  }

  total = 0;
  std::vector<std::pair<double, double>> hops_and_demands;
  for (const auto& hop_and_demand : hop_to_demand) {
    uint32_t hop_count = hop_and_demand.first;
    double demand = hop_and_demand.second;

    total += demand;
    hops_and_demands.emplace_back(hop_count, total);
  }

  cumulative_demands_plot_.AddData("", distances_and_demands);
  cumulative_demands_hop_plot_.AddData("", hops_and_demands);
}

void SingleRCSummaryPlotPack::PlotDemandSizes(
    const nc::lp::DemandMatrix& demand_matrix) {
  const nc::net::GraphStorage* graph = demand_matrix.graph();

  std::vector<std::pair<double, ctr::AggregateId>> decorated;
  for (const auto& element : demand_matrix.elements()) {
    nc::net::GraphNodeIndex src = element.src;
    nc::net::GraphNodeIndex dst = element.dst;
    double demand_Mbps = element.demand.Mbps();
    if (demand_Mbps == 0) {
      continue;
    }

    decorated.push_back(
        std::make_pair(demand_Mbps, ctr::AggregateId(src, dst)));
  }

  std::sort(decorated.begin(), decorated.end());
  ctr::AggregateId max = decorated.back().second;
  ctr::AggregateId med = decorated[decorated.size() / 2].second;

  std::vector<double> demand_sizes;
  for (const auto& value_and_id : decorated) {
    demand_sizes.emplace_back(value_and_id.first);
  }

  demand_sizes_interesting_.Add(
      {"max",
       nc::StrCat(nc::ToStringMaxDecimals(demand_sizes.back(), 2), " Mbps"),
       max.ToString(*graph)});
  demand_sizes_interesting_.Add(
      {"med", nc::StrCat(nc::ToStringMaxDecimals(
                             demand_sizes[demand_sizes.size() / 2], 2),
                         " Mbps"),
       med.ToString(*graph)});
  demand_sizes_plot_.AddData("", demand_sizes);
}

void SingleRCSummaryPlotPack::PlotSPUtilizations(
    const nc::lp::DemandMatrix& demand_matrix) {
  const nc::net::GraphStorage* graph = demand_matrix.graph();

  std::vector<std::pair<double, nc::net::GraphLinkIndex>> decorated;
  nc::net::GraphLinkMap<double> utilizations = demand_matrix.SPUtilization();
  for (const auto& link_and_utilization : utilizations) {
    decorated.emplace_back(*link_and_utilization.second,
                           link_and_utilization.first);
  }

  std::sort(decorated.begin(), decorated.end());
  nc::net::GraphLinkIndex max = decorated.back().second;
  nc::net::GraphLinkIndex med = decorated[decorated.size() / 2].second;

  std::vector<double> link_utilizations;
  for (const auto& value_and_id : decorated) {
    link_utilizations.emplace_back(value_and_id.first);
  }

  sp_utilizations_interesting_.Add(
      {"max", nc::ToStringMaxDecimals(link_utilizations.back(), 2),
       graph->GetLink(max)->ToStringNoPorts()});
  sp_utilizations_interesting_.Add(
      {"med", nc::ToStringMaxDecimals(
                  link_utilizations[link_utilizations.size() / 2], 2),
       graph->GetLink(med)->ToStringNoPorts()});
  sp_utilizations_plot_.AddData("", link_utilizations);
}

void SingleRCSummaryPlotPack::PlotTotalDelayAtLinkScale(
    const nc::lp::DemandMatrix& demand_matrix) {
  const std::map<std::string, std::string>& properties =
      demand_matrix.properties();
  std::string step_str = nc::FindOrDie(properties, "headroom_vs_delay_step");
  std::string delays_at_scale_str =
      nc::FindOrDie(properties, "headroom_vs_delay");

  double abs_stride;
  CHECK(nc::safe_strtod(step_str, &abs_stride));
  const std::vector<double>& delay_at_scale =
      GetCommaSeparated(delays_at_scale_str);

  std::vector<std::pair<double, double>> total_delay_values;
  for (size_t i = 0; i < delay_at_scale.size(); ++i) {
    double multiplier = i * abs_stride;
    total_delay_values.emplace_back(multiplier, delay_at_scale[i]);
  }

  double first_value = total_delay_values.front().second;
  for (auto& multiplier_and_total_delay : total_delay_values) {
    multiplier_and_total_delay.first = 1 - multiplier_and_total_delay.first;
    multiplier_and_total_delay.second =
        multiplier_and_total_delay.second / first_value;
  }
  total_delay_at_link_scale_plot_.AddData("", total_delay_values);
}

void SingleRCSummaryPlotPack::PlotDelayVsUtilization(
    const std::string& opt_label, const RCSummary& rc) {
  nc::net::GraphLinkMap<double> utilizations = rc.utilizations;
  std::vector<std::pair<nc::net::GraphLinkIndex, double>> all_links;

  const nc::net::GraphStorage* graph = rc.parent->demand_matrix->graph();
  for (nc::net::GraphLinkIndex link : graph->AllLinks()) {
    if (utilizations.HasValue(link)) {
      all_links.emplace_back(link, utilizations.GetValueOrDie(link));
    } else {
      all_links.emplace_back(link, 0);
    }
  }

  std::sort(all_links.begin(), all_links.end(),
            [graph](const std::pair<nc::net::GraphLinkIndex, double>& lhs,
                    const std::pair<nc::net::GraphLinkIndex, double>& rhs) {
              nc::net::Delay lhs_delay = graph->GetLink(lhs.first)->delay();
              nc::net::Delay rhs_delay = graph->GetLink(rhs.first)->delay();
              return lhs_delay < rhs_delay;
            });

  std::vector<double> xs;
  std::vector<double> ys;
  for (size_t i = 0; i < all_links.size(); ++i) {
    xs.emplace_back(i);

    double utilization = all_links[i].second;
    ys.emplace_back(utilization);
  }

  link_delay_vs_link_utilization_plot_.AddData(opt_label, xs, ys);
}

void SingleRCSummaryPlotPack::PlotAbsoluteStretch(const std::string& opt_label,
                                                  const RCSummary& rc) {
  nc::DiscreteDistribution<uint64_t> dist;
  for (size_t i = 0; i < rc.abs_stretches.size(); ++i) {
    double rel_stretch = rc.abs_stretches[i];
    double flow_count = rc.flow_counts[i];
    uint64_t value_discrete =
        static_cast<uint64_t>(kDiscreteMultiplier * rel_stretch);
    dist.Add(value_discrete, flow_count);
  }

  std::vector<uint64_t> percentiles = dist.Percentiles(kPercentilesCount);
  CHECK(percentiles.size() == kPercentilesCount + 1);

  std::vector<std::pair<double, double>> path_stretch_abs_data;
  for (size_t i = 0; i < percentiles.size(); ++i) {
    double discrete_value = percentiles[i];
    double p = discrete_value / kDiscreteMultiplier;
    path_stretch_abs_data.emplace_back(
        p, static_cast<double>(i) / kPercentilesCount);
  }

  absolute_path_stretch_plot_.AddData(opt_label, path_stretch_abs_data);
}

}  // namespace alg_eval
}  // namespace ctr
