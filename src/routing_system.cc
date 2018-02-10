#include <gflags/gflags.h>
#include "routing_system.h"

#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "ncode/common.h"
#include "ncode/file.h"
#include "ncode/logging.h"
#include "ncode/lp/demand_matrix.h"
#include "ncode/map_util.h"
#include "ncode/net/net_common.h"
#include "ncode/viz/web_page.h"
#include "metrics/metrics.h"

namespace ctr {

DEFINE_bool(pin_mean, false,
            "If true will always use the mean aggregate level.");

DEFINE_bool(dump_input_to_optimizer, false,
            "If true will dump the inputs to the optimizer.");

static auto* per_aggregate_mean_input =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "aggregate_mean_input",
            "Raw mean level of each aggregate (in Mbps)", "Aggregate source",
            "Aggregate destination");

static auto* per_aggregate_flow_count_input =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string, std::string>(
            "aggregate_flow_count_input",
            "Flow count inputs for each aggregate", "Aggregate source",
            "Aggregate destination");

static auto* per_aggregate_mean_predicted_input =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "aggregate_mean_predicted_input",
            "Predicated mean level of each aggregate (in Mbps)",
            "Aggregate source", "Aggregate destination");

static auto* per_aggregate_predicted_input =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "aggregate_predicted_input",
            "Scaled per-aggregate level. This is the level actually used when "
            "generating the output (in Mbps).",
            "Aggregate source", "Aggregate destination");

static auto* per_aggregate_output_total_delay =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string, std::string>(
            "aggregate_total_delay_micros",
            "Sum of flow_count * path_delay for each of the aggregate's paths "
            "in the output",
            "Aggregate source", "Aggregate destination");

static auto* per_aggregate_output_total_delay_rel =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<uint64_t, std::string, std::string>(
            "aggregate_total_delay_rel_micros",
            "Sum of flow_count * (path_delay - sp_delay) for each of the "
            "aggregate's paths in the output",
            "Aggregate source", "Aggregate destination");

static auto* per_aggregate_output_total_delay_rel_change =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "aggregate_total_delay_rel_change",
            "Sum of flow_count * (path_delay - sp_delay) / sp_delay for each "
            "of the aggregate's paths in the output",
            "Aggregate source", "Aggregate destination");

static auto* per_aggregate_scale_fraction =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "aggregate_scale_fraction",
            "By how much an aggregate had to be scaled", "Aggregate source",
            "Aggregate destination");

static auto* per_link_load =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<double, std::string>(
        "link_load_in_output",
        "What the output of the optimizer says each link should be loaded (in "
        "Mbps)",
        "Link");

static auto* per_link_utilization =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<double, std::string>(
        "link_utilization_in_output",
        "What the output of the optimizer says each link's utilization should "
        "be",
        "Link");

static auto* scale_loop_iteration_count =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<uint64_t>(
        "scale_up_iteration_count",
        "How many iterations of the scale loop happened");

static auto* max_commodity_scale_mean_predicted_input =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<double>(
        "max_commodity_scale_mean_predicted_input",
        "Scale factor for the TM based on mean values after prediction.");

static auto* max_commodity_scale_mean_input =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<double>(
        "max_commodity_scale_mean_input",
        "Scale factor for the TM based on mean values of input.");

static auto* max_commodity_scale_predicted_input =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<double>(
        "max_commodity_scale_predicted_input",
        "Scale factor for the TM used in optimization.");

static auto* path_output_split =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<double, std::string>(
        "path_fraction", "Records output per-path split fractions", "Path");

#define ADD_AGGREGATE_MAP_TO_METRIC(history_map, metric, v)           \
  for (const auto& aggregate_id_and_history : history_map) {          \
    const AggregateId& aggregate_id = aggregate_id_and_history.first; \
    const auto& map_value = aggregate_id_and_history.second;          \
    std::string src_id = graph_->GetNode(aggregate_id.src())->id();   \
    std::string dst_id = graph_->GetNode(aggregate_id.dst())->id();   \
    auto* handle = metric->GetHandle(src_id, dst_id);                 \
    handle->AddValue(v);                                              \
  }

static void UpdateCommodityScaleMetric(
    const std::map<AggregateId, AggregateHistory>& history_map,
    const nc::net::GraphStorage* graph,
    nc::metrics::Metric<double, false>* metric) {
  std::vector<nc::lp::DemandMatrixElement> matrix_elements;
  for (const auto& aggregate_id_and_history : history_map) {
    const AggregateId& aggregate_id = aggregate_id_and_history.first;
    const AggregateHistory& history = aggregate_id_and_history.second;
    matrix_elements.emplace_back(aggregate_id.src(), aggregate_id.dst(),
                                 history.mean_rate());
  }
  nc::lp::DemandMatrix demand_matrix(std::move(matrix_elements), graph);
  metric->GetHandle()->AddValue(demand_matrix.MaxCommodityScaleFactor({}, 1.0));
}

static void RecordPathSplitsAndTotalDelay(
    const RoutingConfiguration& routing_config,
    const nc::net::GraphStorage& graph) {
  nc::net::GraphLinkMap<nc::net::Bandwidth> link_to_load;
  for (nc::net::GraphLinkIndex link : graph.AllLinks()) {
    link_to_load[link] = nc::net::Bandwidth::Zero();
  }

  for (const auto& aggregate_and_routes : routing_config.routes()) {
    const AggregateId& aggregate_id = aggregate_and_routes.first;
    std::string src_id = graph.GetNode(aggregate_id.src())->id();
    std::string dst_id = graph.GetNode(aggregate_id.dst())->id();
    nc::net::Delay sp_delay = aggregate_id.GetSPDelay(graph);

    double total_delay = 0;
    double total_delay_rel = 0;
    double total_delay_rel_change = 0;
    for (const auto& route_and_fraction : aggregate_and_routes.second) {
      const nc::net::Walk* path = route_and_fraction.first;
      nc::net::Delay delay = path->delay();
      double fraction = route_and_fraction.second;

      total_delay += delay.count() * fraction;
      total_delay_rel += (delay.count() - sp_delay.count()) * fraction;
      total_delay_rel_change += ((delay.count() - sp_delay.count()) /
                                 static_cast<double>(sp_delay.count())) *
                                fraction;

      std::string path_string = path->ToStringNoPorts(graph);
      auto* handle = path_output_split->GetHandle(path_string);
      handle->AddValue(fraction);

      const DemandAndFlowCount& initial_demand =
          nc::FindOrDieNoPrint(routing_config.demands(), aggregate_id);
      nc::net::Bandwidth load_on_path = initial_demand.first * fraction;

      for (nc::net::GraphLinkIndex link : path->links()) {
        link_to_load[link] += load_on_path;
      }
    }

    using namespace std::chrono;
    microseconds total_delay_micros = duration_cast<microseconds>(
        nc::net::Delay(static_cast<uint64_t>(total_delay)));
    microseconds total_delay_rel_micros = duration_cast<microseconds>(
        nc::net::Delay(static_cast<uint64_t>(total_delay_rel)));

    per_aggregate_output_total_delay->GetHandle(src_id, dst_id)
        ->AddValue(total_delay_micros.count());
    per_aggregate_output_total_delay_rel->GetHandle(src_id, dst_id)
        ->AddValue(total_delay_rel_micros.count());
    per_aggregate_output_total_delay_rel_change->GetHandle(src_id, dst_id)
        ->AddValue(total_delay_rel_change);
  }

  for (const auto& link_and_load : link_to_load) {
    nc::net::GraphLinkIndex link = link_and_load.first;
    nc::net::Bandwidth load = *link_and_load.second;

    const nc::net::GraphLink* link_ptr = graph.GetLink(link);
    std::string link_id = link_ptr->ToStringNoPorts();
    per_link_load->GetHandle(link_id)->AddValue(load.Mbps());
    per_link_utilization->GetHandle(link_id)
        ->AddValue(load / link_ptr->bandwidth());
  }
}

RoutingSystemUpdateResult RoutingSystem::Update(
    const std::map<AggregateId, AggregateHistory>& history) {
  // First we will predict what the history is expected to be on the next
  // timestep.
  std::map<AggregateId, AggregateHistory> next_history =
      estimator_.EstimateNext(history);
  CHECK(next_history.size() == history.size());

  // The initial input will assign all aggregates to their mean level.
  std::map<AggregateId, DemandAndFlowCount> input =
      GetInitialInput(next_history);

  std::unique_ptr<RoutingConfiguration> output;
  std::map<AggregateId, double> scale_fractions;
  for (const auto& aggregate_and_rest : input) {
    scale_fractions[aggregate_and_rest.first] = 0.0;
  }

  // This will be updated every time the model runs. If the max/mean level are
  // pinned will be empty.
  std::unique_ptr<CompetingAggregates> competing_aggregates;

  size_t i_count = 0;
  while (true) {
    TrafficMatrix tm(graph_, input);

    if (FLAGS_dump_input_to_optimizer) {
      auto demand_matrix = tm.ToDemandMatrix();

      std::string dump_filename =
          nc::Substitute("dump_pass_$0.demands", i_count);
      demand_matrix->ToRepetitaFileOrDie(graph_->NodeOrderOrDie(),
                                         dump_filename);
      LOG(INFO) << "Stored TM for pass " << i_count << " at " << dump_filename;
    }

    output = optimizer_->Optimize(tm);
    CHECK(output->demands().size() == tm.demands().size());
    ++i_count;

    // After optimizing we should check to see which aggregates go over links
    // that do not fit.
    std::set<AggregateId> aggregates_no_fit;
    std::tie(aggregates_no_fit, competing_aggregates) =
        CheckWithProbModel(*output, next_history);
    LOG(INFO) << aggregates_no_fit.size() << " aggregates do not fit at pass "
              << i_count;
    if (aggregates_no_fit.empty()) {
      break;
    }

    if (FLAGS_pin_mean) {
      competing_aggregates->Clear();
      break;
    }

    // Need to scale up the aggregates that do not fit.
    bool scaled_any = ScaleUpAggregates(aggregates_no_fit, next_history, &input,
                                        &scale_fractions);
    if (!scaled_any) {
      // All aggregates at their max.
      break;
    }
  }

  CHECK(output);
  if (config_.store_to_metrics) {
    // Will first record the raw input.
    ADD_AGGREGATE_MAP_TO_METRIC(history, per_aggregate_mean_input,
                                map_value.mean_rate().Mbps());
    UpdateCommodityScaleMetric(history, graph_, max_commodity_scale_mean_input);

    // And flow counts.
    ADD_AGGREGATE_MAP_TO_METRIC(history, per_aggregate_flow_count_input,
                                map_value.flow_count());

    UpdateCommodityScaleMetric(next_history, graph_,
                               max_commodity_scale_mean_predicted_input);

    // Record the predicted mean level.
    ADD_AGGREGATE_MAP_TO_METRIC(next_history,
                                per_aggregate_mean_predicted_input,
                                map_value.mean_rate().Mbps());

    // Record the number of iterations.
    scale_loop_iteration_count->GetHandle()->AddValue(i_count);

    double scale_factor =
        output->ToDemandMatrix()->MaxCommodityScaleFactor({}, 1.0);
    max_commodity_scale_predicted_input->GetHandle()->AddValue(scale_factor);

    // And also the per-path splits from the output.
    RecordPathSplitsAndTotalDelay(*output, *graph_);

    // Record by how much we had to scale each aggregate to make it fit.
    ADD_AGGREGATE_MAP_TO_METRIC(scale_fractions, per_aggregate_scale_fraction,
                                map_value);

    // Will record the input that generated the output.
    ADD_AGGREGATE_MAP_TO_METRIC(input, per_aggregate_predicted_input,
                                map_value.first.Mbps());
  }

  return {std::move(output), std::move(competing_aggregates), next_history};
}

std::pair<std::set<AggregateId>, std::unique_ptr<CompetingAggregates>>
RoutingSystem::CheckWithProbModel(
    const RoutingConfiguration& routing,
    const std::map<AggregateId, AggregateHistory>& histories) {
  CHECK(histories.size() == routing.routes().size());

  nc::net::GraphLinkMap<std::vector<std::pair<AggregateId, double>>>
      per_link_aggregates;
  for (const auto& aggregate_and_routes : routing.routes()) {
    const AggregateId& aggregate = aggregate_and_routes.first;
    const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;
    for (const auto& route_and_fraction : routes) {
      const nc::net::Walk* walk = route_and_fraction.first;
      double fraction = route_and_fraction.second;

      for (nc::net::GraphLinkIndex link : walk->links()) {
        per_link_aggregates[link].emplace_back(aggregate, fraction);
      }
    }
  }

  std::vector<ProbModelQuery> queries;
  std::vector<nc::net::GraphLinkIndex> links_in_order;
  for (const auto& link_and_aggregates : per_link_aggregates) {
    nc::net::GraphLinkIndex link = link_and_aggregates.first;
    const std::vector<std::pair<AggregateId, double>>&
        aggregates_and_fractions = *link_and_aggregates.second;
    links_in_order.emplace_back(link);

    queries.emplace_back();
    ProbModelQuery& query = queries.back();
    query.type = ProbModelQuery::BOTH;
    query.rate = graph_->GetLink(link)->bandwidth();
    query.aggregates = aggregates_and_fractions;
  }

  ProbModel prob_model(config_.prob_model_config);
  for (const auto& aggregate_and_history : histories) {
    const AggregateId& aggregate = aggregate_and_history.first;
    const AggregateHistory& history = aggregate_and_history.second;
    prob_model.AddAggregate(aggregate, &history);
  }
  std::vector<ProbModelReply> replies = prob_model.Query(queries);

  // Will need the replies indexed by link, so that we can later figure out
  // which link along a path is most likely to congest.
  nc::net::GraphLinkMap<ProbModelReply> replies_by_link;

  std::set<AggregateId> aggregates_no_fit;
  for (size_t i = 0; i < replies.size(); ++i) {
    nc::net::GraphLinkIndex link = links_in_order[i];
    const ProbModelReply& reply = replies[i];
    replies_by_link[link] = reply;
    if (reply.fit) {
      continue;
    }

    const std::vector<std::pair<AggregateId, double>>&
        aggregates_and_fractions = per_link_aggregates.GetValueOrDie(link);
    for (const auto& aggregate_and_fraction : aggregates_and_fractions) {
      aggregates_no_fit.insert(aggregate_and_fraction.first);
    }
  }

  auto competing_aggregates = nc::make_unique<CompetingAggregates>();
  for (const auto& aggregate_and_routes : routing.routes()) {
    const AggregateId& aggregate = aggregate_and_routes.first;

    // We only care about aggregates that can fit the demand. If there are one
    // or many links along any of the aggregate's paths that cannot fit the
    // demand, we will add no state for the aggregate.
    bool aggregate_fits = true;

    // State to add to competing_aggregates for this aggregate.
    std::vector<AggregatesAndCapacity> to_add;

    const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;
    for (const auto& route_and_fraction : routes) {
      const nc::net::Walk* walk = route_and_fraction.first;

      double max_capacity_fraction = -1;
      nc::net::GraphLinkIndex link_most_likely_to_congest;

      // Need to find the most likely to congest link along the path, and
      // add the aggregates that cross it as competing aggregates to the current
      // aggregate.
      for (nc::net::GraphLinkIndex link : walk->links()) {
        const ProbModelReply& reply = replies_by_link.GetValueOrDie(link);
        if (!reply.fit) {
          aggregate_fits = false;
          break;
        }

        double fraction;
        if (reply.optimal_rate == nc::net::Bandwidth::Zero()) {
          // The path fits the traffic trivially.
          fraction = 0;
        } else {
          nc::net::Bandwidth link_capacity = graph_->GetLink(link)->bandwidth();
          fraction = reply.optimal_rate / link_capacity;
        }

        if (fraction > max_capacity_fraction) {
          link_most_likely_to_congest = link;
          max_capacity_fraction = fraction;
        }
      }

      if (!aggregate_fits) {
        break;
      }

      CHECK(max_capacity_fraction != -1);
      to_add.emplace_back();
      AggregatesAndCapacity& aggregates_and_capacity = to_add.back();
      aggregates_and_capacity.capacity =
          graph_->GetLink(link_most_likely_to_congest)->bandwidth();
      for (const auto& aggregate_and_fraction :
           per_link_aggregates.GetValueOrDie(link_most_likely_to_congest)) {
        // In competing_states we do not want to add aggregates that compete
        // with themselves.
        if (aggregate_and_fraction.first == aggregate) {
          continue;
        }

        aggregates_and_capacity.aggregates.emplace_back(aggregate_and_fraction);
      }
    }

    if (!aggregate_fits) {
      continue;
    }

    competing_aggregates->AddAggregatesAndCapacity(aggregate, to_add);
  }

  return {aggregates_no_fit, std::move(competing_aggregates)};
}

bool RoutingSystem::ScaleUpAggregates(
    const std::set<AggregateId>& aggregates,
    const std::map<AggregateId, AggregateHistory>& histories,
    std::map<AggregateId, DemandAndFlowCount>* out,
    std::map<AggregateId, double>* scale_fractions) {
  CHECK(!aggregates.empty());
  bool scaled_at_least_one = false;
  for (const auto& aggregate_and_history : histories) {
    const AggregateId& aggregate_id = aggregate_and_history.first;
    if (!nc::ContainsKey(aggregates, aggregate_id)) {
      continue;
    }

    const AggregateHistory& history = aggregate_and_history.second;
    nc::net::Bandwidth mean_rate = history.mean_rate();
    nc::net::Bandwidth max_rate = history.max_rate();

    const double range_mbps = max_rate.Mbps() - mean_rate.Mbps();
    if (range_mbps == 0) {
      continue;
    }

    DemandAndFlowCount& current_demand_and_flow_count =
        nc::FindOrDieNoPrint(*out, aggregate_id);
    nc::net::Bandwidth& current_rate = current_demand_and_flow_count.first;
    CHECK(current_rate >= mean_rate);

    const double mean_rate_mbps = mean_rate.Mbps();
    const double current_rate_mbps = current_rate.Mbps();
    const double curr_scale = current_rate_mbps / mean_rate_mbps;

    // Will scale the aggregate's rate to be somewhere between the mean and
    // the max. Need to first figure out where it is now.
    double curr_fraction = mean_rate_mbps * (curr_scale - 1) / range_mbps;
    if (curr_fraction > 0.99) {
      continue;
    }

    // Want to avoid updating all_fit in the case where all aggregates at
    // already at 1.0.
    scaled_at_least_one = true;

    // Will move each aggregate 10% closer to its max rate.
    double new_fraction = curr_fraction + kScaleFraction;
    new_fraction = std::min(new_fraction, 1.0);
    (*scale_fractions)[aggregate_id] = new_fraction;

    double new_rate_mbps = mean_rate_mbps + range_mbps * new_fraction;
    current_rate = nc::net::Bandwidth::FromMBitsPerSecond(new_rate_mbps);
  }

  return scaled_at_least_one;
}

std::map<AggregateId, DemandAndFlowCount> RoutingSystem::GetInitialInput(
    const std::map<AggregateId, AggregateHistory>& histories) {
  std::map<AggregateId, DemandAndFlowCount> out;
  for (const auto& aggregate_and_history : histories) {
    const AggregateId& aggregate_id = aggregate_and_history.first;
    const AggregateHistory& history = aggregate_and_history.second;

    nc::net::Bandwidth init_rate = history.mean_rate();
    init_rate =
        std::max(init_rate, nc::net::Bandwidth::FromKBitsPerSecond(10.0));
    out[aggregate_id] = {init_rate, history.flow_count()};
  }

  return out;
}

}  // namespace ctr
