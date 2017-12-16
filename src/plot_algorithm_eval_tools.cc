#include "plot_algorithm_eval_tools.h"

#include <algorithm>
#include <functional>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/stats.h"
#include "metrics/metrics_parser.h"

static constexpr char kAggregateSPDelayMetric[] = "aggregate_sp_delay_ms";
static constexpr char kAggregateRateMetric[] = "aggregate_rate_Mbps";
static constexpr char kAggregatePathCountMetric[] = "opt_path_count";
static constexpr char kPathDelayMetric[] = "opt_path_delay_ms";
static constexpr char kPathCountMetric[] = "opt_path_flow_count";
static constexpr char kLinkUtilizationMetric[] = "opt_link_utilization";
static constexpr size_t kDiscreteMultiplier = 10000;
static constexpr size_t kPercentilesCount = 10000;

namespace ctr {
namespace alg_eval {

std::vector<AggregateTMState> OptimizerTMState::GetAggregates() const {
  size_t path_index = 0;
  std::vector<AggregateTMState> out;
  for (size_t i = 0; i < aggregate_path_count.size(); ++i) {
    uint32_t path_count = aggregate_path_count[i];
    uint32_t sp_delay_ms = parent_state->aggregate_sp_delay_ms[i];
    double rate_Mbps = parent_state->aggregate_rate_Mbps[i];

    std::vector<std::pair<uint32_t, uint32_t>> paths;
    for (size_t j = 0; j < path_count; ++j) {
      uint32_t delay_ms = path_delay_ms[path_index];
      uint32_t flow_count = path_flow_count[path_index];
      paths.emplace_back(delay_ms, flow_count);

      ++path_index;
    }

    out.emplace_back();
    out.back().paths = std::move(paths);
    out.back().sp_delay_ms = sp_delay_ms;
    out.back().rate_Mbps = rate_Mbps;
  }

  return out;
}

OptimizerTMState* DataStorage::GetOptimizerTMState(const std::string& topology,
                                                   const std::string& optimizer,
                                                   const std::string& tm) {
  std::unique_lock<std::mutex> lock(mu_);
  TMStateMap& tm_state_map = data_[optimizer];
  std::unique_ptr<OptimizerTMState>& opt_tm_state_ptr =
      tm_state_map[{topology, tm}];

  if (!opt_tm_state_ptr) {
    opt_tm_state_ptr = nc::make_unique<OptimizerTMState>();
    opt_tm_state_ptr->parent_state = GetTMStatePrivate(topology, tm);
  }
  return opt_tm_state_ptr.get();
}

TMState* DataStorage::GetTMState(const std::string& topology,
                                 const std::string& tm) {
  std::unique_lock<std::mutex> lock(mu_);
  return GetTMStatePrivate(topology, tm);
}

TMState* DataStorage::GetTMStatePrivate(const std::string& topology,
                                        const std::string& tm) {
  std::unique_ptr<TMState>& tm_state_ptr = tm_state_data_[{topology, tm}];
  if (!tm_state_ptr) {
    tm_state_ptr = nc::make_unique<TMState>();
  }
  return tm_state_ptr.get();
}

class SingleMetricProcessor : public nc::metrics::parser::MetricProcessor {
 public:
  using UpdateStateCallback = std::function<void(
      const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state)>;

  SingleMetricProcessor(const std::string& metric, DataStorage* storage,
                        const std::set<TopologyAndTM>* to_ignore,
                        UpdateStateCallback callback);

  bool InterestedInMetric(const std::string& metric_id) override {
    return metric_ == metric_id;
  }

  bool InterestedInFields(const nc::metrics::PBManifestEntry& manifest_entry,
                          uint32_t manifest_index) override;

  void ProcessEntry(const nc::metrics::PBMetricEntry& entry,
                    const nc::metrics::PBManifestEntry& manifest_entry,
                    uint32_t manifest_index) override;

 private:
  const std::string metric_;

  UpdateStateCallback callback_;

  // Combinations of TM and topology to ignore.
  const std::set<TopologyAndTM>* to_ignore_;

  // The storage that produces OptimizerTMState instances.
  DataStorage* storage_;

  // Maps from a manifest index to state.
  std::map<uint32_t, OptimizerTMState*> manifest_index_to_state_;
};

class SingleTMMetricProcessor : public nc::metrics::parser::MetricProcessor {
 public:
  using UpdateStateCallback = std::function<void(
      const nc::metrics::PBMetricEntry& entry, TMState* tm_state)>;

  SingleTMMetricProcessor(const std::string& metric, DataStorage* storage,
                          const std::set<TopologyAndTM>* to_ignore,
                          UpdateStateCallback callback);

  bool InterestedInMetric(const std::string& metric_id) override {
    return metric_ == metric_id;
  }

  bool InterestedInFields(const nc::metrics::PBManifestEntry& manifest_entry,
                          uint32_t manifest_index) override;

  void ProcessEntry(const nc::metrics::PBMetricEntry& entry,
                    const nc::metrics::PBManifestEntry& manifest_entry,
                    uint32_t manifest_index) override {
    nc::Unused(manifest_entry);
    TMState* tm_state = nc::FindOrDie(manifest_index_to_state_, manifest_index);
    callback_(entry, tm_state);
  }

 private:
  const std::string metric_;

  UpdateStateCallback callback_;

  // Combinations of TM and topology to ignore.
  const std::set<TopologyAndTM>* to_ignore_;

  // The storage that produces OptimizerTMState instances.
  DataStorage* storage_;

  // Maps from a manifest index to state.
  std::map<uint32_t, TMState*> manifest_index_to_state_;
};

SingleMetricProcessor::SingleMetricProcessor(
    const std::string& metric, DataStorage* storage,
    const std::set<TopologyAndTM>* to_ignore, UpdateStateCallback callback)
    : metric_(metric),
      callback_(callback),
      to_ignore_(to_ignore),
      storage_(storage) {}

bool SingleMetricProcessor::InterestedInFields(
    const nc::metrics::PBManifestEntry& manifest_entry,
    uint32_t manifest_index) {
  // Fields are assumed to be topology,TM,optimizer.
  CHECK(manifest_entry.fields_size() == 3) << metric_;
  CHECK(manifest_entry.fields(0).type() ==
        nc::metrics::PBMetricField_Type_STRING);
  CHECK(manifest_entry.fields(1).type() ==
        nc::metrics::PBMetricField_Type_STRING);
  CHECK(manifest_entry.fields(2).type() ==
        nc::metrics::PBMetricField_Type_STRING);

  const std::string& topology = manifest_entry.fields(0).string_value();
  const std::string& tm = manifest_entry.fields(1).string_value();
  const std::string& optimizer = manifest_entry.fields(2).string_value();

  if (nc::ContainsKey(*to_ignore_, make_pair(topology, tm))) {
    return false;
  }

  OptimizerTMState* tm_state =
      storage_->GetOptimizerTMState(topology, optimizer, tm);

  CHECK(!nc::ContainsKey(manifest_index_to_state_, manifest_index));
  manifest_index_to_state_[manifest_index] = tm_state;
  return true;
}

void SingleMetricProcessor::ProcessEntry(
    const nc::metrics::PBMetricEntry& entry,
    const nc::metrics::PBManifestEntry& manifest_entry,
    uint32_t manifest_index) {
  nc::Unused(manifest_entry);
  OptimizerTMState* tm_state =
      nc::FindOrDie(manifest_index_to_state_, manifest_index);
  callback_(entry, tm_state);
}

SingleTMMetricProcessor::SingleTMMetricProcessor(
    const std::string& metric, DataStorage* storage,
    const std::set<TopologyAndTM>* to_ignore, UpdateStateCallback callback)
    : metric_(metric),
      callback_(callback),
      to_ignore_(to_ignore),
      storage_(storage) {}

bool SingleTMMetricProcessor::InterestedInFields(
    const nc::metrics::PBManifestEntry& manifest_entry,
    uint32_t manifest_index) {
  // Fields are assumed to be topology,TM.
  CHECK(manifest_entry.fields_size() == 2) << metric_;
  CHECK(manifest_entry.fields(0).type() ==
        nc::metrics::PBMetricField_Type_STRING);
  CHECK(manifest_entry.fields(1).type() ==
        nc::metrics::PBMetricField_Type_STRING);

  const std::string& topology = manifest_entry.fields(0).string_value();
  const std::string& tm = manifest_entry.fields(1).string_value();
  if (nc::ContainsKey(*to_ignore_, make_pair(topology, tm))) {
    return false;
  }

  TMState* tm_state = storage_->GetTMState(topology, tm);

  CHECK(!nc::ContainsKey(manifest_index_to_state_, manifest_index));
  manifest_index_to_state_[manifest_index] = tm_state;
  return true;
}

std::pair<std::vector<double>, std::vector<double>> GetStretchDistribution(
    const TMStateMap& tm_state_map) {
  nc::DiscreteDistribution<uint64_t> dist;
  std::vector<double> maxs;

  for (const auto& key_and_data : tm_state_map) {
    const OptimizerTMState& tm_state = *(key_and_data.second);
    std::vector<AggregateTMState> aggregates = tm_state.GetAggregates();
    double max = 0;
    for (const auto& aggregate : aggregates) {
      double sp_delay_ms = aggregate.sp_delay_ms;
      for (const auto& delay_and_count : aggregate.paths) {
        double path_delay_ms = delay_and_count.first;
        uint32_t flow_count = delay_and_count.second;
        CHECK(path_delay_ms >= sp_delay_ms);
        double abs_stretch = path_delay_ms - sp_delay_ms;
        double rel_stretch = abs_stretch / sp_delay_ms;
        max = std::max(rel_stretch, max);

        // Have to discretize the value.
        uint64_t value_discrete =
            static_cast<uint64_t>(kDiscreteMultiplier * abs_stretch);
        dist.Add(value_discrete, flow_count);
      }
    }
    maxs.emplace_back(max);
  }

  std::vector<uint64_t> percentiles = dist.Percentiles(kPercentilesCount);
  CHECK(percentiles.size() == kPercentilesCount + 1);

  std::vector<double> out;
  out.reserve(percentiles.size());
  for (uint64_t discrete_value : percentiles) {
    out.emplace_back(static_cast<double>(discrete_value) / kDiscreteMultiplier);
  }

  return {out, maxs};
}

std::vector<double> GetStretchDistribution(const OptimizerTMState& tm_state,
                                           bool absolute) {
  nc::DiscreteDistribution<uint64_t> dist;
  std::vector<AggregateTMState> aggregates = tm_state.GetAggregates();
  for (const auto& aggregate : aggregates) {
    double sp_delay_ms = aggregate.sp_delay_ms;
    for (const auto& delay_and_count : aggregate.paths) {
      double path_delay_ms = delay_and_count.first;
      uint32_t flow_count = delay_and_count.second;
      CHECK(path_delay_ms >= sp_delay_ms);
      double stretch;
      if (absolute) {
        stretch = path_delay_ms - sp_delay_ms;
      } else {
        stretch = (path_delay_ms - sp_delay_ms) / sp_delay_ms;
      }

      // Have to discretize the value.
      uint64_t value_discrete =
          static_cast<uint64_t>(kDiscreteMultiplier * stretch);
      dist.Add(value_discrete, flow_count);
    }
  }

  std::vector<uint64_t> percentiles = dist.Percentiles(kPercentilesCount);
  CHECK(percentiles.size() == kPercentilesCount + 1);

  std::vector<double> out;
  out.reserve(percentiles.size());
  for (uint64_t discrete_value : percentiles) {
    out.emplace_back(static_cast<double>(discrete_value) / kDiscreteMultiplier);
  }

  return out;
}

std::unique_ptr<DataStorage> ParseTMGenUtilMetrics(
    const std::string& location, const std::set<TopologyAndTM>& to_ignore) {
  auto data_storage = nc::make_unique<DataStorage>();
  auto sp_delay_processor = nc::make_unique<SingleTMMetricProcessor>(
      kAggregateSPDelayMetric, data_storage.get(), &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, TMState* tm_state) {
        tm_state->aggregate_sp_delay_ms.emplace_back(entry.uint32_value());
      });

  auto rate_processor = nc::make_unique<SingleTMMetricProcessor>(
      kAggregateRateMetric, data_storage.get(), &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, TMState* tm_state) {
        tm_state->aggregate_rate_Mbps.emplace_back(entry.double_value());
      });

  auto path_delay_processor = nc::make_unique<SingleMetricProcessor>(
      kPathDelayMetric, data_storage.get(), &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->path_delay_ms.emplace_back(entry.uint32_value());
      });

  auto path_flow_count_processor = nc::make_unique<SingleMetricProcessor>(
      kPathCountMetric, data_storage.get(), &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->path_flow_count.emplace_back(entry.uint32_value());
      });

  auto link_utilization_processor = nc::make_unique<SingleMetricProcessor>(
      kLinkUtilizationMetric, data_storage.get(), &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->link_utilization.emplace_back(entry.double_value());
      });

  auto aggregate_path_count_processor = nc::make_unique<SingleMetricProcessor>(
      kAggregatePathCountMetric, data_storage.get(), &to_ignore,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->aggregate_path_count.emplace_back(entry.uint32_value());
      });

  nc::metrics::parser::MetricsParser parser(location);
  parser.AddProcessor(std::move(sp_delay_processor));
  parser.AddProcessor(std::move(rate_processor));
  parser.AddProcessor(std::move(path_delay_processor));
  parser.AddProcessor(std::move(path_flow_count_processor));
  parser.AddProcessor(std::move(link_utilization_processor));
  parser.AddProcessor(std::move(aggregate_path_count_processor));
  parser.Parse();

  return data_storage;
}

}  // namespace alg_eval
}  // namespace ctr
