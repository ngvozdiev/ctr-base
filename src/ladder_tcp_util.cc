// A ladder-like topology with three aggregates.

#include <gflags/gflags.h>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/htsim/match.h"
#include "ncode_common/src/htsim/tcp.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "common.h"
#include "controller.h"
#include "mean_est/mean_est.h"
#include "metrics/metrics.h"
#include "net_instrument.h"
#include "net_mock.h"
#include "opt/ctr.h"
#include "opt/opt.h"
#include "opt/path_provider.h"
#include "routing_system.h"
#include "tldr.h"

using namespace std::chrono;
using namespace ctr::controller;
using namespace nc::net;

DEFINE_uint64(tcp_short_flow_object_size_bytes, 50000ul,
              "Size of each object that the short flows will transmit.");
DEFINE_uint64(
    tcp_short_flow_wait_time_ms, 5000,
    "Average time a short-flow source will wait before transmission.");
DEFINE_uint64(tcp_short_flow_access_link_rate_Mbps, 10,
              "Access link rate for short flows");
DEFINE_uint64(tcp_short_flow_count, 10, "Number of short-flow clients");
DEFINE_bool(only_short_flows, false,
            "If true only short TCP flows will be added");
DEFINE_uint64(increase_start_ms, 60 * 1000 * 5,
              "When to start increasing traffic");
DEFINE_uint64(increase_end_ms, 60 * 1000 * 10,
              "When to stop increasing traffic");
DEFINE_uint64(decrease_start_ms, 60 * 1000 * 20,
              "When to start decreasing traffic");
DEFINE_uint64(decrease_end_ms, 60 * 1000 * 25,
              "When to stop decreasing traffic");
DEFINE_string(opt, "CTR", "The optimizer to use");
DEFINE_double(decay_factor, 0.0, "How quickly to decay prediction");
DEFINE_uint64(period_duration_ms, 60000, "Length of the period");
DEFINE_uint64(history_bin_size_ms, 100, "How big each history bin is");
DEFINE_bool(simulate_initial_handshake, true,
            "Whether or not to simulate an initial handshake");
DEFINE_uint64(
    tcp_initial_cwnd_pkts, 4,
    "How many packets should there be in the inital congestion window");
DEFINE_uint64(duration_sec, 500, "For how long to run (in simulated time)");

static void AddKeyFrames(Bandwidth start_rate, Bandwidth mid_rate,
                         Bandwidth end_rate, size_t count,
                         TCPFlowGroup* tcp_flow_group) {
  using namespace std::chrono;
  tcp_flow_group->AddKeyFrame(milliseconds(FLAGS_increase_start_ms), count,
                              start_rate);
  tcp_flow_group->AddKeyFrame(milliseconds(FLAGS_increase_end_ms), count,
                              mid_rate);
  tcp_flow_group->AddKeyFrame(milliseconds(FLAGS_decrease_start_ms), count,
                              mid_rate);
  tcp_flow_group->AddKeyFrame(milliseconds(FLAGS_decrease_end_ms), count,
                              end_rate);
}

static void SetConstantRate(Bandwidth rate, size_t count,
                            TCPFlowGroup* tcp_flow_group) {
  tcp_flow_group->AddKeyFrame(milliseconds(60 * 1000 * 5), count, rate);
}

static std::pair<TCPFlowGroup, TCPFlowGroup> AddFlows(size_t num_flows) {
  using namespace std::chrono;
  TCPFlowGroup short_flow_group(
      FLAGS_tcp_short_flow_count, milliseconds(5), milliseconds(20),
      FLAGS_tcp_short_flow_object_size_bytes,
      milliseconds(FLAGS_tcp_short_flow_wait_time_ms));
  short_flow_group.set_mean_object_size_fixed(true);
  short_flow_group.set_mean_wait_time_fixed(false);
  short_flow_group.set_access_link_rate_spread(0.3);
  short_flow_group.set_random_access_link_queue(false);
  short_flow_group.set_initial_time_offset(milliseconds(10000));

  short_flow_group.AddKeyFrame(
      milliseconds(10), FLAGS_tcp_short_flow_count,
      Bandwidth::FromMBitsPerSecond(FLAGS_tcp_short_flow_access_link_rate_Mbps *
                                    FLAGS_tcp_short_flow_count));

  TCPFlowGroup long_flow_group(num_flows, milliseconds(5), milliseconds(20),
                               std::numeric_limits<uint64_t>::max(),
                               milliseconds::zero());
  long_flow_group.set_mean_object_size_fixed(true);
  long_flow_group.set_mean_wait_time_fixed(true);
  long_flow_group.set_access_link_rate_spread(0.3);
  long_flow_group.set_random_access_link_queue(true);
  return {short_flow_group, long_flow_group};
}

static constexpr size_t kLongFlowCount = 500;

static void AddTopStream(NetworkContainer* container) {
  GraphNodeIndex src = container->graph()->NodeFromStringOrDie("N0");
  GraphNodeIndex dst = container->graph()->NodeFromStringOrDie("N1");
  container->AddAggregate(
      {src, dst}, ctr::GetDummyHistory(Bandwidth::FromMBitsPerSecond(500),
                                       milliseconds(FLAGS_history_bin_size_ms),
                                       milliseconds(FLAGS_period_duration_ms),
                                       kLongFlowCount));
  container->AddAggregate(
      {dst, src}, ctr::GetDummyHistory(Bandwidth::FromMBitsPerSecond(5),
                                       milliseconds(FLAGS_history_bin_size_ms),
                                       milliseconds(FLAGS_period_duration_ms),
                                       kLongFlowCount));

  // Stream between N0 and N1. Initially it will be small, but it will grow
  // later. When it grows it should "push" down the other aggregates.
  auto short_and_long_flow_groups = AddFlows(kLongFlowCount);
  TCPFlowGroup& short_flow_group = short_and_long_flow_groups.first;
  TCPFlowGroup& long_flow_group = short_and_long_flow_groups.second;

  if (!FLAGS_only_short_flows) {
    // Flows will start at 500mbps, double and then shrink back.
    AddKeyFrames(
        Bandwidth::FromMBitsPerSecond(500), Bandwidth::FromMBitsPerSecond(1250),
        Bandwidth::FromMBitsPerSecond(500), kLongFlowCount, &long_flow_group);
  } else {
    AddKeyFrames(
        Bandwidth::FromMBitsPerSecond(5), Bandwidth::FromMBitsPerSecond(12),
        Bandwidth::FromMBitsPerSecond(5), kLongFlowCount, &long_flow_group);
  }

  container->AddTCPFlowGroup({src, dst}, short_flow_group);
  container->AddTCPFlowGroup({src, dst}, long_flow_group);
}

static void AddBottomStream(const std::string& src_id,
                            const std::string& dst_id,
                            NetworkContainer* container) {
  GraphNodeIndex src = container->graph()->NodeFromStringOrDie(src_id);
  GraphNodeIndex dst = container->graph()->NodeFromStringOrDie(dst_id);

  // Flows will be constantly at 800Mbps.
  Bandwidth total_rate = FLAGS_only_short_flows
                             ? Bandwidth::FromMBitsPerSecond(8)
                             : Bandwidth::FromMBitsPerSecond(800);

  container->AddAggregate(
      {src, dst}, ctr::GetDummyHistory(
                      total_rate, milliseconds(FLAGS_history_bin_size_ms),
                      milliseconds(FLAGS_period_duration_ms), kLongFlowCount));

  // Need another aggregate for ACKs.
  container->AddAggregate(
      {dst, src}, ctr::GetDummyHistory(total_rate / 100.0,
                                       milliseconds(FLAGS_history_bin_size_ms),
                                       milliseconds(FLAGS_period_duration_ms),
                                       kLongFlowCount));

  auto short_and_long_flow_groups = AddFlows(kLongFlowCount);
  TCPFlowGroup& short_flow_group = short_and_long_flow_groups.first;
  TCPFlowGroup& long_flow_group = short_and_long_flow_groups.second;
  SetConstantRate(total_rate, kLongFlowCount, &long_flow_group);
  container->AddTCPFlowGroup({src, dst}, short_flow_group);
  container->AddTCPFlowGroup({src, dst}, long_flow_group);
}

// A global variable that will keep a reference to the event queue, useful for
// logging, as the logging handler only accepts a C-style function pointer.
static nc::EventQueue* event_queue_global_ptr;

void CustomLogHandler(nc::LogLevel level, const char* filename, int line,
                      const std::string& message, nc::LogColor color) {
  uint64_t time_as_ms = event_queue_global_ptr->TimeToRawMillis(
      event_queue_global_ptr->CurrentTime());
  std::string new_message = nc::StrCat(time_as_ms, "ms ", message);
  nc::DefaultLogHandler(level, filename, line, new_message, color);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::metrics::InitMetrics();

  // The event queue.
  nc::SimTimeEventQueue event_queue;
  event_queue_global_ptr = &event_queue;
  nc::SetLogHandler(CustomLogHandler);

  auto timestamp_provider =
      ::nc::make_unique<nc::metrics::SimTimestampProvider>(&event_queue);
  nc::metrics::DefaultMetricManager()->set_timestamp_provider(
      std::move(timestamp_provider));

  nc::net::GraphBuilder builder = nc::net::GenerateLadder(
      4, Bandwidth::FromGBitsPerSecond(10.0), milliseconds(5), 0.1,
      {milliseconds(1), milliseconds(5), milliseconds(20), milliseconds(50)});
  nc::net::GraphStorage graph(builder);

  ctr::PathProvider path_provider(&graph);
  std::unique_ptr<ctr::Optimizer> opt;
  if (FLAGS_opt == "CTR") {
    opt = nc::make_unique<ctr::CTROptimizer>(&path_provider, 0.95, true, false);
  } else if (FLAGS_opt == "B4") {
    opt = nc::make_unique<ctr::B4Optimizer>(&path_provider, false, 0.95);
  } else if (FLAGS_opt == "B4(P)") {
    opt = nc::make_unique<ctr::B4Optimizer>(&path_provider, true, 0.95);
  } else if (FLAGS_opt == "MinMax") {
    opt = nc::make_unique<ctr::MinMaxOptimizer>(&path_provider, 0.9, false);
  } else if (FLAGS_opt == "SP") {
    opt = nc::make_unique<ctr::ShortestPathOptimizer>(&path_provider);
  }

  ctr::MeanScaleEstimatorFactory estimator_factory(
      {1.05, FLAGS_decay_factor, FLAGS_decay_factor, 10});
  ctr::RoutingSystemConfig routing_system_config;
  ctr::RoutingSystem routing_system(routing_system_config, opt.get(),
                                    &estimator_factory);

  nc::net::IPAddress controller_ip(100);
  nc::net::IPAddress device_ip_base(20000);
  nc::net::IPAddress tldr_ip_base(30000);
  nc::net::IPAddress flow_group_ip_base(40000);
  nc::net::DevicePortNumber enter_port(4000);
  nc::net::DevicePortNumber exit_port(5000);

  std::chrono::milliseconds tcp_rto_timer_period(10);
  std::chrono::milliseconds round_duration(FLAGS_period_duration_ms);
  std::chrono::milliseconds poll_period(FLAGS_history_bin_size_ms);

  ctr::controller::Controller controller(controller_ip, &routing_system,
                                         &event_queue, &graph);
  ctr::controller::NetworkContainerConfig containter_config(
      device_ip_base, tldr_ip_base, flow_group_ip_base, enter_port, exit_port,
      std::chrono::milliseconds(100), tcp_rto_timer_period, false);
  containter_config.tcp_config.inital_cwnd_size = FLAGS_tcp_initial_cwnd_pkts;
  containter_config.tcp_config.simulate_initial_handshake =
      FLAGS_simulate_initial_handshake;

  nc::ThresholdEnforcerPolicy te_policy;
  ctr::TLDRConfig tldr_config(te_policy, {}, nc::htsim::kWildIPAddress,
                              nc::htsim::kWildIPAddress, controller_ip,
                              round_duration, poll_period, 100);
  ctr::controller::NetworkContainer network_container(
      containter_config, tldr_config, &graph, &controller, &event_queue);

  // Records per-path stats.
  ctr::InputPacketObserver packet_observer(network_container.controller(),
                                           std::chrono::milliseconds(10),
                                           &event_queue);

  ctr::DefaultDeviceFactory device_factory;
  network_container.AddElementsFromGraph(&device_factory, &packet_observer);

  AddTopStream(&network_container);
  AddBottomStream("N2", "N3", &network_container);
  AddBottomStream("N6", "N7", &network_container);
  AddBottomStream("N10", "N11", &network_container);

  network_container.InitAggregatesInController();

  ctr::NetInstrument net_instrument(network_container.internal_queues(),
                                    network_container.flow_group_tcp_sources(),
                                    std::chrono::milliseconds(10),
                                    &event_queue);
  nc::htsim::ProgressIndicator progress_indicator(
      std::chrono::milliseconds(100), &event_queue);
  event_queue.RunAndStopIn(std::chrono::seconds(FLAGS_duration_sec));
  nc::metrics::DefaultMetricManager()->PersistAllMetrics();
}
