#include "routing_system.h"

#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/net_common.h"
#include "metrics/metrics.h"

namespace ctr {

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

#define ADD_AGGREGATE_MAP_TO_METRIC(history_map, metric, v)             \
  if (metric_timestamp != std::numeric_limits<uint64_t>::max()) {       \
    for (const auto& aggregate_id_and_history : history_map) {          \
      const AggregateId& aggregate_id = aggregate_id_and_history.first; \
      const auto& map_value = aggregate_id_and_history.second;          \
      std::string src_id = graph_->GetNode(aggregate_id.src())->id();   \
      std::string dst_id = graph_->GetNode(aggregate_id.dst())->id();   \
      auto* handle = metric->GetHandle(src_id, dst_id);                 \
      handle->AddValueWithTimestamp(metric_timestamp, v);               \
    }                                                                   \
  }

static void UpdateCommodityScaleMetric(
    uint64_t time_now,
    const std::map<AggregateId, AggregateHistory>& history_map,
    const nc::net::GraphStorage* graph,
    nc::metrics::Metric<double, false>* metric) {
  if (time_now == std::numeric_limits<uint64_t>::max()) {
    return;
  }

  std::vector<nc::lp::DemandMatrixElement> matrix_elements;
  for (const auto& aggregate_id_and_history : history_map) {
    const AggregateId& aggregate_id = aggregate_id_and_history.first;
    const AggregateHistory& history = aggregate_id_and_history.second;
    matrix_elements.emplace_back(aggregate_id.src(), aggregate_id.dst(),
                                 history.mean_rate());
  }
  nc::lp::DemandMatrix demand_matrix(std::move(matrix_elements), graph);
  metric->GetHandle()->AddValueWithTimestamp(
      time_now, demand_matrix.MaxCommodityScaleFractor());
}

static void RecordPathSplitsAndTotalDelay(
    const RoutingConfiguration& routing_config,
    const nc::net::GraphStorage& graph, uint64_t timestamp) {
  if (timestamp == std::numeric_limits<uint64_t>::max()) {
    return;
  }

  for (const auto& aggregate_and_routes : routing_config.routes()) {
    const AggregateId& aggregate_id = aggregate_and_routes.first;
    std::string src_id = graph.GetNode(aggregate_id.src())->id();
    std::string dst_id = graph.GetNode(aggregate_id.dst())->id();

    nc::net::Delay sp_delay = aggregate_id.GetSPDelay(graph);
    const DemandAndFlowCount& demand_and_flow_count =
        nc::FindOrDieNoPrint(routing_config.demands(), aggregate_id);
    size_t flow_count = demand_and_flow_count.second;

    double total_delay = 0;
    double total_delay_rel = 0;
    double total_delay_rel_change = 0;
    for (const auto& route_and_fraction : aggregate_and_routes.second) {
      const nc::net::Walk* path = route_and_fraction.first;
      nc::net::Delay delay = path->delay();

      total_delay += delay.count() * flow_count;
      total_delay_rel += (delay.count() - sp_delay.count()) * flow_count;
      total_delay_rel_change += ((delay.count() - sp_delay.count()) /
                                 static_cast<double>(sp_delay.count())) *
                                flow_count;

      std::string path_string = path->ToStringNoPorts(graph);
      double fraction = route_and_fraction.second;
      auto* handle = path_output_split->GetHandle(path_string);
      handle->AddValueWithTimestamp(timestamp, fraction);
    }

    using namespace std::chrono;
    microseconds total_delay_micros = duration_cast<microseconds>(
        nc::net::Delay(static_cast<uint64_t>(total_delay)));
    microseconds total_delay_rel_micros = duration_cast<microseconds>(
        nc::net::Delay(static_cast<uint64_t>(total_delay_rel)));

    per_aggregate_output_total_delay->GetHandle(src_id, dst_id)
        ->AddValueWithTimestamp(timestamp, total_delay_micros.count());
    per_aggregate_output_total_delay_rel->GetHandle(src_id, dst_id)
        ->AddValueWithTimestamp(timestamp, total_delay_rel_micros.count());
    per_aggregate_output_total_delay_rel_change->GetHandle(src_id, dst_id)
        ->AddValueWithTimestamp(timestamp, total_delay_rel_change);
  }
}

std::unique_ptr<RoutingConfiguration> RoutingSystem::Update(
    const std::map<AggregateId, AggregateHistory>& history,
    uint64_t metric_timestamp) {
  // Will first record the raw input.
  ADD_AGGREGATE_MAP_TO_METRIC(history, per_aggregate_mean_input,
                              map_value.mean_rate().Mbps());
  UpdateCommodityScaleMetric(metric_timestamp, history, graph_,
                             max_commodity_scale_mean_input);

  // And flow counts.
  ADD_AGGREGATE_MAP_TO_METRIC(history, per_aggregate_flow_count_input,
                              map_value.flow_count());

  // First we will predict what the history is expected to be on the next
  // timestep.
  std::map<AggregateId, AggregateHistory> next_history =
      estimator_.EstimateNext(history);
  UpdateCommodityScaleMetric(metric_timestamp, next_history, graph_,
                             max_commodity_scale_mean_predicted_input);

  // Record the predicted mean level.
  ADD_AGGREGATE_MAP_TO_METRIC(next_history, per_aggregate_mean_predicted_input,
                              map_value.mean_rate().Mbps());

  // The initial input will assign all aggregates to their mean level.
  std::map<AggregateId, DemandAndFlowCount> input =
      GetInitialInput(next_history);

  std::unique_ptr<RoutingConfiguration> output;
  std::map<AggregateId, double> scale_fractions;
  for (const auto& aggregate_and_rest : input) {
    scale_fractions[aggregate_and_rest.first] = 0.0;
  }

  size_t i_count = 0;
  while (true) {
    TrafficMatrix tm(graph_, input);
    output = optimizer_->Optimize(tm);
    ++i_count;

    // After optimizing we should check to see which aggregates go over links
    // that do not fit.
    std::set<AggregateId> aggregates_no_fit =
        CheckWithProbModel(*output, next_history);
    if (aggregates_no_fit.empty()) {
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
  if (metric_timestamp != std::numeric_limits<uint64_t>::max()) {
    // Record the number of iterations.
    scale_loop_iteration_count->GetHandle()->AddValueWithTimestamp(
        metric_timestamp, i_count);

    double scale_factor = output->ToDemandMatrix()->MaxCommodityScaleFractor();
    max_commodity_scale_predicted_input->GetHandle()->AddValueWithTimestamp(
        metric_timestamp, scale_factor);

    // And also the per-path splits from the output.
    RecordPathSplitsAndTotalDelay(*output, *graph_, metric_timestamp);
  }

  // Record by how much we had to scale each aggregate to make it fit.
  ADD_AGGREGATE_MAP_TO_METRIC(scale_fractions, per_aggregate_scale_fraction,
                              map_value);

  // Will record the input that generated the output.
  ADD_AGGREGATE_MAP_TO_METRIC(input, per_aggregate_predicted_input,
                              map_value.first.Mbps());

  return output;
}

std::set<AggregateId> RoutingSystem::CheckWithProbModel(
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
  std::vector<const std::vector<std::pair<AggregateId, double>>*> map_values;
  for (const auto& link_and_aggregates : per_link_aggregates) {
    nc::net::GraphLinkIndex link = link_and_aggregates.first;
    const std::vector<std::pair<AggregateId, double>>&
        aggregates_and_fractions = *link_and_aggregates.second;
    map_values.emplace_back(&aggregates_and_fractions);

    queries.emplace_back();
    ProbModelQuery& query = queries.back();
    query.type = ProbModelQuery::BOTH;
    query.rate = graph_->GetLink(link)->bandwidth();
    query.aggregates = aggregates_and_fractions;
  }

  ProbModel prob_model(prob_model_config_);
  for (const auto& aggregate_and_history : histories) {
    const AggregateId& aggregate = aggregate_and_history.first;
    const AggregateHistory& history = aggregate_and_history.second;
    prob_model.AddAggregate(aggregate, &history);
  }
  std::vector<ProbModelReply> replies = prob_model.Query(queries);

  std::set<AggregateId> out;
  for (size_t i = 0; i < replies.size(); ++i) {
    const ProbModelReply& reply = replies[i];
    if (reply.fit) {
      continue;
    }

    for (const auto& aggregate_and_fraction : *(map_values[i])) {
      out.insert(aggregate_and_fraction.first);
    }
  }

  return out;
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

    out[aggregate_id] = {history.mean_rate(), history.flow_count()};
  }

  return out;
}

}  // namespace ctr
