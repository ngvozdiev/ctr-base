#include "net_mock.h"

#include <algorithm>
#include <limits>
#include <utility>

#include "ncode_common/src/map_util.h"
#include "routing_system.h"
#include "metrics/metrics.h"

namespace ctr {

static auto* link_utilization_metric =
    nc::metrics::DefaultMetricManager()
        -> GetUnsafeMetric<double, std::string, std::string>(
            "link_utilization", "Records per-link utilization", "Link source",
            "Link destination");

static auto* path_split_metric =
    nc::metrics::DefaultMetricManager() -> GetUnsafeMetric<double, std::string>(
        "path_fraction", "Records per-path split fractions", "Path");

NetMock::NetMock(const std::map<AggregateId, BinSequence>& initial_sequences,
                 std::chrono::milliseconds period_duration,
                 std::chrono::milliseconds history_bin_size,
                 RoutingSystem* routing_system)
    : history_bin_size_(history_bin_size),
      initial_sequences_(initial_sequences),
      routing_system_(routing_system),
      graph_(routing_system_->graph()) {
  CHECK(!initial_sequences.empty());
  size_t min_bin_count = std::numeric_limits<size_t>::max();
  std::chrono::milliseconds bin_size = std::chrono::milliseconds::zero();
  for (const auto& id_and_bin_sequence : initial_sequences) {
    const BinSequence& bin_sequence = id_and_bin_sequence.second;
    if (bin_size == std::chrono::milliseconds::zero()) {
      bin_size = std::chrono::duration_cast<std::chrono::milliseconds>(
          bin_sequence.bin_size());
    } else {
      CHECK(bin_size == bin_sequence.bin_size());
    }

    size_t count = bin_sequence.bin_count();
    min_bin_count = std::min(min_bin_count, count);
  }

  period_duration_bins_ = period_duration.count() / bin_size.count();
  CHECK(period_duration_bins_ > 0);
  period_count_ = min_bin_count / period_duration_bins_;
}

// Generates the input to the system.
std::map<AggregateId, AggregateHistory> NetMock::GenerateInput(
    std::map<AggregateId, BinSequence>* period_sequences) const {
  std::map<AggregateId, AggregateHistory> input;
  for (auto& aggregate_and_bins : *period_sequences) {
    const AggregateId& aggregate = aggregate_and_bins.first;
    BinSequence& bins = aggregate_and_bins.second;

    input.emplace(
        std::piecewise_construct, std::forward_as_tuple(aggregate),
        std::forward_as_tuple(bins.GenerateHistory(history_bin_size_)));
  }

  return input;
}

nc::net::GraphLinkMap<std::vector<double>> NetMock::CheckOutput(
    const std::map<AggregateId, BinSequence>& period_sequences,
    const RoutingConfiguration& configuration) const {
  nc::net::GraphLinkMap<std::unique_ptr<BinSequence>> link_to_bins;

  // First need to figure out which paths cross each link. Will also build a
  // map from paths to aggregates and path indices.
  for (const auto& aggregate_and_routes : configuration.routes()) {
    const AggregateId& aggregate = aggregate_and_routes.first;
    const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;

    std::vector<double> fractions;
    for (const auto& route_and_fraction : routes) {
      fractions.emplace_back(route_and_fraction.second);
    }

    // For each of the aggregate's paths, the bins that go on that path.
    std::vector<BinSequence> aggregate_split =
        nc::FindOrDieNoPrint(period_sequences, aggregate).SplitOrDie(fractions);

    for (size_t i = 0; i < routes.size(); ++i) {
      const nc::net::Walk* path = routes[i].first;
      for (nc::net::GraphLinkIndex link : path->links()) {
        std::unique_ptr<BinSequence>& bin_sequence_ptr = link_to_bins[link];
        if (!bin_sequence_ptr) {
          bin_sequence_ptr = nc::make_unique<BinSequence>(aggregate_split[i]);
        } else {
          bin_sequence_ptr->Combine(aggregate_split[i]);
        }
      }
    }
  }

  nc::net::GraphLinkMap<std::vector<double>> out;
  for (const auto& link_and_bins : link_to_bins) {
    nc::net::GraphLinkIndex link = link_and_bins.first;
    nc::net::Bandwidth rate = graph_->GetLink(link)->bandwidth();
    out[link] = (*link_and_bins.second)->Residuals(rate);
  }

  return out;
}

std::map<AggregateId, BinSequence> NetMock::GetNthPeriod(size_t n) const {
  size_t period_start_bin = n * period_duration_bins_;
  size_t period_end_bin = (n + 1) * period_duration_bins_;

  std::map<AggregateId, BinSequence> out;
  for (const auto& aggregate_and_bins : initial_sequences_) {
    const AggregateId& aggregate = aggregate_and_bins.first;
    const BinSequence& bins = aggregate_and_bins.second;

    BinSequence to_end = bins.CutFromStart(period_end_bin);
    BinSequence period_sequence = to_end.Offset(period_start_bin);
    out.emplace(std::piecewise_construct, std::forward_as_tuple(aggregate),
                std::forward_as_tuple(period_sequence));
  }

  return out;
}

std::unique_ptr<RoutingConfiguration> NetMock::InitialOutput() const {
  std::map<AggregateId, BinSequence> zero_period = GetNthPeriod(0);
  std::map<AggregateId, AggregateHistory> input = GenerateInput(&zero_period);
  return routing_system_->Update(input);
}

static size_t CheckSameSize(
    const nc::net::GraphLinkMap<std::vector<double>>& values) {
  size_t i = 0;
  for (const auto& link_and_values : values) {
    const std::vector<double>& v = *link_and_values.second;
    CHECK(v.size() > 0);
    if (i == 0) {
      i = v.size();
    } else {
      CHECK(i == v.size());
    }
  }

  return i;
}

void NetMock::RecordPathSplits(const RoutingConfiguration& routing_config,
                             uint64_t timestamp) const {
  for (const auto& aggregate_and_routes : routing_config.routes()) {
    for (const auto& route_and_fraction : aggregate_and_routes.second) {
      std::string path = route_and_fraction.first->ToStringNoPorts(*graph_);
      double fraction = route_and_fraction.second;
      auto* handle = path_split_metric->GetHandle(path);
      handle->AddValueWithTimestamp(timestamp, fraction);
    }
  }
}

void NetMock::Run() {
  std::unique_ptr<RoutingConfiguration> output = InitialOutput();
  RecordPathSplits(*output, 0);
  size_t timestamp = 0;
  for (size_t i = 0; i < period_count_; ++i) {
    LOG(ERROR) << "Period " << i;

    std::map<AggregateId, BinSequence> period_sequences = GetNthPeriod(i);
    nc::net::GraphLinkMap<std::vector<double>> per_link_residuals =
        CheckOutput(period_sequences, *output);
    size_t num_residuals = CheckSameSize(per_link_residuals);

    for (auto link_and_residuals : per_link_residuals) {
      nc::net::GraphLinkIndex link = link_and_residuals.first;
      std::vector<double>& residuals = *link_and_residuals.second;

      const nc::net::GraphLink* link_ptr = graph_->GetLink(link);
      auto* handle = link_utilization_metric->GetHandle(link_ptr->src_id(),
                                                        link_ptr->dst_id());
      size_t t = timestamp;
      for (double v : residuals) {
        handle->AddValueWithTimestamp(t++, v);
      }
    }

    timestamp += num_residuals;
    std::map<AggregateId, AggregateHistory> input =
        GenerateInput(&period_sequences);
    output = routing_system_->Update(input);
    RecordPathSplits(*output, timestamp);
  }
}

}  // namespace ctr
