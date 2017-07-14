#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/viz/grapher.h"
#include "../common.h"
#include "../mean_est/mean_est.h"
#include "../routing_system.h"
#include "ctr.h"
#include "path_provider.h"

DEFINE_double(cd_link_gbps, 1.0, "Capacity of C->D");
DEFINE_double(ab_link_gbps, 1.0, "Capacity of A->B");
DEFINE_double(cd_link_ms, 10, "Delay of C->D");
DEFINE_double(ab_link_ms, 10, "Delay of A->B");
DEFINE_double(cd_aggregate_gbps, 0.9, "Demand of C->D aggregate");
DEFINE_double(ab_aggregate_gbps, 0.9, "Demand of A->B aggregate");
DEFINE_uint64(cd_aggregate_flows, 1000, "Number of C->D flows");
DEFINE_uint64(ab_aggregate_flows, 1000, "Number of A->B flows");
DEFINE_uint64(steps, 20, "Number of steps");
DEFINE_double(decay_factor, 0.0, "How quickly to decay prediction");
DEFINE_string(opt, "CTR", "The optimizer to use");
DEFINE_bool(dampen_ratio, false, "Preserve the flow count / volume ratio?");

namespace ctr {

static constexpr std::chrono::milliseconds kSyntheticBinSize =
    std::chrono::milliseconds(10);
static constexpr size_t kSyntheticBinCount = 600;

class StabilityEvalHarness {
 public:
  StabilityEvalHarness(TrafficMatrix* initial_tm, RoutingSystem* routing_system,
                       size_t cycle_count)
      : initial_tm_(initial_tm),
        routing_system_(routing_system),
        graph_(routing_system_->graph()) {
    Cycle(cycle_count);
  }

  const std::vector<RoutingConfigurationDelta>& deltas() const {
    return deltas_;
  }

  void PlotLinkFractions(const std::string& out) const {
    std::vector<nc::viz::DataSeries2D> series;
    for (const auto& path_and_fractions : path_to_fraction_) {
      const nc::net::Walk* walk = path_and_fractions.first;

      series.emplace_back();
      nc::viz::DataSeries2D& data_series = series.back();
      data_series.data = path_and_fractions.second;
      data_series.label = walk->ToStringNoPorts(*graph_);
    }

    nc::viz::PlotParameters2D params;
    params.x_label = "timestep";
    params.y_label = "fraction";

    nc::viz::PythonGrapher grapher(out);
    grapher.PlotLine(params, series);
  }

  void DumpLinkFractions() const {
    std::cout << path_to_fraction_.size() << "\n";
    for (const auto& path_and_fractions : path_to_fraction_) {
      const nc::net::Walk* walk = path_and_fractions.first;

      std::vector<std::string> to_combine;
      for (const auto& times_and_fractions : path_and_fractions.second) {
        to_combine.emplace_back(nc::StrCat("(", times_and_fractions.first, ",",
                                           times_and_fractions.second, ")"));
      }

      std::cout << walk->ToStringNoPorts(*graph_) << " : ("
                << nc::Join(to_combine, ",") << ")\n";
    }
  }

  void DumpAggregateVolumes() const {
    std::cout << volumes_.size() << "\n";
    for (const auto& aggregate_and_volumes : volumes_) {
      const AggregateId& aggregate_id = aggregate_and_volumes.first;
      const std::vector<nc::net::Bandwidth>& volumes =
          aggregate_and_volumes.second;

      std::string out = nc::Join(volumes, ",", [](nc::net::Bandwidth v) {
        return std::to_string(v.Mbps());
      });
      std::cout << aggregate_id.ToString(*graph_) << " : (" << out << ")\n";
    }
  }

 private:
  // Produces a traffic matrix that is scaled based on each flow's shortest path
  // stretch. If all flows of an aggregate are on a path that is 2x as long as
  // the shortest path then their demand will be 1/2 of what it would be if they
  // were on the shortest path.
  std::unique_ptr<TrafficMatrix> ScaleBasedOnOutput(
      const RoutingConfiguration& routing,
      const TrafficMatrix& original_tm) const {
    std::map<AggregateId, DemandAndFlowCount> out;
    for (const auto& aggregate_and_routes : routing.routes()) {
      const AggregateId& aggregate_id = aggregate_and_routes.first;

      const DemandAndFlowCount& demand_and_flow_count =
          nc::FindOrDieNoPrint(routing.demands(), aggregate_id);
      const DemandAndFlowCount& initial_demand_and_flow_count =
          nc::FindOrDieNoPrint(initial_tm_->demands(), aggregate_id);

      size_t flow_count = demand_and_flow_count.second;
      nc::net::Bandwidth sp_bandwidth = initial_demand_and_flow_count.first;
      nc::net::Delay sp_delay = aggregate_id.GetSPDelay(*graph_);

      double total_mbps = 0;
      double per_flow_sp_mbps = sp_bandwidth.Mbps() / flow_count;
      const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;
      for (const RouteAndFraction& route_and_fraction : routes) {
        nc::net::Delay path_delay = route_and_fraction.first->delay();
        double sp_ratio =
            static_cast<double>(sp_delay.count()) / path_delay.count();

        // Each flow in the aggregate will get bandwidth proportional to how far
        // away it is from the shortest path.
        total_mbps += route_and_fraction.second * flow_count * sp_ratio *
                      per_flow_sp_mbps;
        CHECK(total_mbps == total_mbps);
      }

      std::mt19937 rnd(1);
      std::uniform_real_distribution<double> dist(1, 10);

      double new_flow_count = flow_count;
      if (FLAGS_dampen_ratio) {
        // The demand in the output will not be the same as the one in the input
        // to the optimizer because of prediction. Need to dampen the ratio
        // based on the original input, not the predicted one.
        const DemandAndFlowCount& original_demand_and_flow_count =
            nc::FindOrDieNoPrint(original_tm.demands(), aggregate_id);

        new_flow_count = total_mbps * flow_count /
                         original_demand_and_flow_count.first.Mbps();
        new_flow_count = std::max(1.0, new_flow_count);
      }

      out[aggregate_id] = {nc::net::Bandwidth::FromMBitsPerSecond(total_mbps),
                           new_flow_count};
    }

    return nc::make_unique<TrafficMatrix>(graph_, out);
  }

  // Performs a number of runs with the optimizer. Returns n - 1 values, one for
  // each delta between a run and its following run.
  void Cycle(size_t n) {
    std::unique_ptr<TrafficMatrix> prev_tm;
    std::unique_ptr<RoutingConfiguration> prev_routing;
    for (size_t i = 0; i < n; ++i) {
      const TrafficMatrix* tm = prev_tm ? prev_tm.get() : initial_tm_;
      UpdateVolumes(*tm);

      std::unique_ptr<RoutingConfiguration> routing;
      std::tie(routing, std::ignore) =
          routing_system_->Update(SyntheticHistoryFromTM(*tm));

      for (const auto& aggregate_and_routes : routing->routes()) {
        for (const auto& path_and_fraction : aggregate_and_routes.second) {
          path_to_fraction_[path_and_fraction.first].emplace_back(
              i, path_and_fraction.second);
        }
      }

      if (prev_routing) {
        deltas_.emplace_back(prev_routing->GetDifference(*routing));
      }
      prev_tm = ScaleBasedOnOutput(*routing, *tm);
      prev_routing = std::move(routing);
    }
  }

  void UpdateVolumes(const TrafficMatrix& tm) {
    for (const auto& aggregate_and_demands : tm.demands()) {
      const AggregateId& aggregate = aggregate_and_demands.first;
      volumes_[aggregate].emplace_back(aggregate_and_demands.second.first);
    }
  }

  std::map<AggregateId, AggregateHistory> SyntheticHistoryFromTM(
      const TrafficMatrix& tm) {
    std::map<AggregateId, AggregateHistory> out;
    for (const auto& aggregate_demand_and_flow_count : tm.demands()) {
      const AggregateId& aggregate = aggregate_demand_and_flow_count.first;
      nc::net::Bandwidth demand = aggregate_demand_and_flow_count.second.first;
      size_t flow_count = aggregate_demand_and_flow_count.second.second;

      out.emplace(std::piecewise_construct, std::forward_as_tuple(aggregate),
                  std::forward_as_tuple(demand, kSyntheticBinCount,
                                        kSyntheticBinSize, flow_count));
    }

    return out;
  }

  // In the initial TM each aggregate's demand is what it would be if it were
  // routed on its shortest path.
  TrafficMatrix* initial_tm_;

  // The optimizer.
  RoutingSystem* routing_system_;

  // The graph.
  const nc::net::GraphStorage* graph_;

  // N - 1 deltas for each step to the next one.
  std::vector<RoutingConfigurationDelta> deltas_;

  // Per-aggregate volumes.
  std::map<AggregateId, std::vector<nc::net::Bandwidth>> volumes_;

  // Path to fractions of capacity.
  std::map<const nc::net::Walk*, std::vector<std::pair<double, double>>>
      path_to_fraction_;
};

static void RunWithSimpleTopologyTwoAggregates() {
  nc::net::GraphBuilder builder;
  std::chrono::microseconds ab_delay(
      static_cast<size_t>(FLAGS_ab_link_ms * 1000));
  std::chrono::microseconds cd_delay(
      static_cast<size_t>(FLAGS_cd_link_ms * 1000));

  builder.AddLink({"A", "B",
                   nc::net::Bandwidth::FromGBitsPerSecond(FLAGS_ab_link_gbps),
                   ab_delay});
  builder.AddLink({"B", "A",
                   nc::net::Bandwidth::FromGBitsPerSecond(FLAGS_ab_link_gbps),
                   ab_delay});
  builder.AddLink({"A", "C", nc::net::Bandwidth::FromGBitsPerSecond(1),
                   std::chrono::milliseconds(1)});
  builder.AddLink({"C", "A", nc::net::Bandwidth::FromGBitsPerSecond(1),
                   std::chrono::milliseconds(1)});
  builder.AddLink({"C", "D",
                   nc::net::Bandwidth::FromGBitsPerSecond(FLAGS_cd_link_gbps),
                   cd_delay});
  builder.AddLink({"D", "C",
                   nc::net::Bandwidth::FromGBitsPerSecond(FLAGS_cd_link_gbps),
                   cd_delay});
  builder.AddLink({"B", "D", nc::net::Bandwidth::FromGBitsPerSecond(1),
                   std::chrono::milliseconds(1)});
  builder.AddLink({"D", "B", nc::net::Bandwidth::FromGBitsPerSecond(1),
                   std::chrono::milliseconds(1)});

  nc::net::GraphStorage graph(builder);
  TrafficMatrix initial_tm(&graph);
  AggregateId id_one(
      {graph.NodeFromStringOrDie("A"), graph.NodeFromStringOrDie("B")});
  AggregateId id_two(
      {graph.NodeFromStringOrDie("C"), graph.NodeFromStringOrDie("D")});

  initial_tm.AddDemand(
      id_one, {nc::net::Bandwidth::FromGBitsPerSecond(FLAGS_ab_aggregate_gbps),
               FLAGS_ab_aggregate_flows});
  initial_tm.AddDemand(
      id_two, {nc::net::Bandwidth::FromGBitsPerSecond(FLAGS_cd_aggregate_gbps),
               FLAGS_cd_aggregate_flows});

  PathProvider path_provider(&graph);
  std::unique_ptr<Optimizer> opt;
  if (FLAGS_opt == "CTR") {
    opt = nc::make_unique<CTROptimizer>(&path_provider, 1.0, true);
  } else if (FLAGS_opt == "B4") {
    opt = nc::make_unique<B4Optimizer>(&path_provider, false, 1.0);
  } else if (FLAGS_opt == "B4(P)") {
    opt = nc::make_unique<B4Optimizer>(&path_provider, true, 1.0);
  } else if (FLAGS_opt == "MinMax") {
    opt = nc::make_unique<MinMaxOptimizer>(&path_provider, 1.0);
  }

  MeanScaleEstimatorFactory estimator_factory(
      {1.1, FLAGS_decay_factor, FLAGS_decay_factor, 10});
  RoutingSystemConfig routing_system_config;
  routing_system_config.store_to_metrics = false;
  RoutingSystem routing_system(routing_system_config, opt.get(),
                               &estimator_factory);
  StabilityEvalHarness harness(&initial_tm, &routing_system, FLAGS_steps);
  harness.DumpAggregateVolumes();
  harness.DumpLinkFractions();
}

}  // namespace ctr

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  //  ctr::RunWithSimpleTopology();
  ctr::RunWithSimpleTopologyTwoAggregates();

  return 0;
}
