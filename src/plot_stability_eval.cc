#include <gflags/gflags.h>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <mutex>

#include "ncode/common.h"
#include "ncode/logging.h"
#include "ncode/map_util.h"
#include "ncode/strutil.h"
#include "ncode/viz/grapher.h"
#include "../build/metrics.pb.h"
#include "metrics/metrics_parser.h"

DEFINE_string(metrics_dir, "", "The metrics directory.");
DEFINE_string(optimizers, "MinMaxLD,CTR_LIM,B4,CTR,MinMaxK10",
              "Optimizers to plot.");

static constexpr char kVolumeChangeMetric[] = "stability_volume_delta";
static constexpr char kVolumeChangeLongerPathMetric[] =
    "stability_volume_delta_longer_path";
static constexpr char kAddCountMetric[] = "stability_add_count";
static constexpr char kRemoveCountMetric[] = "stability_remove_count";
static constexpr char kUpdateCountMetric[] = "stability_update_count";
static constexpr char kDelayDeltaMetric[] = "stability_total_flow_delay_delta";

using namespace nc::metrics::parser;

// State for a single traffic matrix.
struct OptimizerTMState {
  double volume_delta = 0;
  double volume_on_longer_path = 0;
  double delay_delta = 0;
  size_t add_count = 0;
  size_t update_count = 0;
  size_t remove_count = 0;
};

using TopologyAndTM = std::pair<std::string, std::string>;
using TMStateMap = std::map<TopologyAndTM, std::unique_ptr<OptimizerTMState>>;

using DataVector = std::vector<double>;
using DataMap = std::map<std::pair<std::string, std::string>, DataVector>;

class DataStorage {
 public:
  OptimizerTMState* GetTMState(const std::string& topology,
                               const std::string& optimizer,
                               const std::string& tm) {
    std::unique_lock<std::mutex> lock(mu_);
    TMStateMap& tm_state_map = data_[optimizer];
    std::unique_ptr<OptimizerTMState>& tm_state_ptr =
        tm_state_map[{topology, tm}];

    if (!tm_state_ptr) {
      tm_state_ptr = nc::make_unique<OptimizerTMState>();
    }
    return tm_state_ptr.get();
  }

  const std::map<std::string, TMStateMap>& data() const { return data_; }

 private:
  std::map<std::string, TMStateMap> data_;

  // Protects 'data_'.
  std::mutex mu_;
};

class SingleMetricProcessor : public MetricProcessor {
 public:
  using UpdateStateCallback = std::function<void(
      const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state)>;

  SingleMetricProcessor(const std::string& metric, DataStorage* storage,
                        UpdateStateCallback callback)
      : metric_(metric), callback_(callback), storage_(storage) {}

  bool InterestedInMetric(const std::string& metric_id) override {
    return metric_ == metric_id;
  }

  bool InterestedInFields(const nc::metrics::PBManifestEntry& manifest_entry,
                          uint32_t manifest_index) override {
    // Fields are assumed to be topology,TM,optimizer.
    CHECK(manifest_entry.fields_size() == 3);
    CHECK(manifest_entry.fields(0).type() ==
          nc::metrics::PBMetricField_Type_STRING);
    CHECK(manifest_entry.fields(1).type() ==
          nc::metrics::PBMetricField_Type_STRING);
    CHECK(manifest_entry.fields(2).type() ==
          nc::metrics::PBMetricField_Type_STRING);

    const std::string& topology = manifest_entry.fields(0).string_value();
    const std::string& tm = manifest_entry.fields(1).string_value();
    const std::string& optimizer = manifest_entry.fields(2).string_value();

    OptimizerTMState* tm_state = storage_->GetTMState(topology, optimizer, tm);

    CHECK(!nc::ContainsKey(manifest_index_to_state_, manifest_index));
    manifest_index_to_state_[manifest_index] = tm_state;
    return true;
  }

  void ProcessEntry(const nc::metrics::PBMetricEntry& entry,
                    const nc::metrics::PBManifestEntry& manifest_entry,
                    uint32_t manifest_index) override {
    nc::Unused(manifest_entry);
    OptimizerTMState* tm_state =
        nc::FindOrDie(manifest_index_to_state_, manifest_index);
    callback_(entry, tm_state);
  }

 private:
  const std::string metric_;

  UpdateStateCallback callback_;

  // The storage that produces OptimizerTMState instances.
  DataStorage* storage_;

  // Maps from a manifest index to state.
  std::map<uint32_t, OptimizerTMState*> manifest_index_to_state_;
};

static void PlotCDF(const std::vector<nc::viz::DataSeries1D>& to_plot,
                    const std::string& output) {
  nc::viz::CDFPlot plot;
  plot.AddData(to_plot);
  plot.PlotToDir(output);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_metrics_dir.empty()) << "need --metrics_dir";

  DataStorage data_storage;
  auto volume_processor = nc::make_unique<SingleMetricProcessor>(
      kVolumeChangeMetric, &data_storage,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->volume_delta = entry.double_value();
      });

  auto volume_longer_path_processor = nc::make_unique<SingleMetricProcessor>(
      kVolumeChangeLongerPathMetric, &data_storage,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->volume_on_longer_path = entry.double_value();
      });

  auto add_count_processor = nc::make_unique<SingleMetricProcessor>(
      kAddCountMetric, &data_storage,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->add_count = entry.uint32_value();
      });

  auto remove_count_processor = nc::make_unique<SingleMetricProcessor>(
      kRemoveCountMetric, &data_storage,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->remove_count = entry.uint32_value();
      });

  auto update_count_processor = nc::make_unique<SingleMetricProcessor>(
      kUpdateCountMetric, &data_storage,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->update_count = entry.uint32_value();
      });

  auto delay_delta_processor = nc::make_unique<SingleMetricProcessor>(
      kDelayDeltaMetric, &data_storage,
      [](const nc::metrics::PBMetricEntry& entry, OptimizerTMState* tm_state) {
        tm_state->delay_delta = entry.double_value();
      });

  MetricsParser parser(FLAGS_metrics_dir);
  parser.AddProcessor(std::move(volume_processor));
  parser.AddProcessor(std::move(volume_longer_path_processor));
  parser.AddProcessor(std::move(add_count_processor));
  parser.AddProcessor(std::move(remove_count_processor));
  parser.AddProcessor(std::move(update_count_processor));
  parser.AddProcessor(std::move(delay_delta_processor));
  parser.Parse();

  std::vector<std::string> optimizers = nc::Split(FLAGS_optimizers, ",");
  const std::map<std::string, TMStateMap>& data = data_storage.data();

  std::vector<nc::viz::DataSeries1D> volume_to_plot;
  std::vector<nc::viz::DataSeries1D> volume_on_shorter_path_to_plot;
  std::vector<nc::viz::DataSeries1D> num_paths_updated_to_plot;
  nc::viz::DataSeries1D delay_delta_to_plot;

  for (const std::string& opt : optimizers) {
    LOG(INFO) << "Handle " << opt;
    volume_to_plot.push_back({opt, {}});
    volume_on_shorter_path_to_plot.push_back({opt, {}});
    num_paths_updated_to_plot.push_back({opt, {}});

    nc::viz::DataSeries1D& volume_data_series = volume_to_plot.back();
    nc::viz::DataSeries1D& volume_on_shorter_path_data_series =
        volume_on_shorter_path_to_plot.back();
    nc::viz::DataSeries1D& num_paths_updated_data_series =
        num_paths_updated_to_plot.back();

    const TMStateMap& state_map = nc::FindOrDie(data, opt);
    for (const auto& topolog_and_tm_and_rest : state_map) {
      const OptimizerTMState& tm_state = *(topolog_and_tm_and_rest.second);
      volume_data_series.data.emplace_back(tm_state.volume_delta);
      volume_on_shorter_path_data_series.data.emplace_back(
          tm_state.volume_delta - tm_state.volume_on_longer_path);
      num_paths_updated_data_series.data.emplace_back(
          tm_state.add_count + tm_state.remove_count + tm_state.update_count);

      if (opt == "CTR_LIM") {
        delay_delta_to_plot.data.emplace_back(tm_state.delay_delta);
      }
    }
  }

  PlotCDF(volume_to_plot, "stability_total_volume_out");
  PlotCDF(volume_on_shorter_path_to_plot,
          "stability_total_volume_shorter_path_out");
  PlotCDF(num_paths_updated_to_plot, "stability_num_paths_updated_out");
  PlotCDF({delay_delta_to_plot}, "stability_delay_delta_out");
}
