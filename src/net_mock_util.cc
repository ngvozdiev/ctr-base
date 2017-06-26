#include <gflags/gflags.h>
#include <chrono>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "common.h"
#include "mean_est/mean_est.h"
#include "net_mock.h"
#include "opt/opt.h"
#include "opt/ctr.h"
#include "pcap_data.h"
#include "routing_system.h"
#include "metrics/metrics.h"

DEFINE_string(topology, "", "A file with a topology");
DEFINE_string(traffic_matrix, "", "A file with a traffic matrix");
DEFINE_string(pcap_trace_store, "trace_store.pb",
              "A file with information about .pcap traces");
DEFINE_uint64(period_duration_ms, 60000, "Length of the period");
DEFINE_uint64(history_bin_size_ms, 100, "How big each history bin is");
DEFINE_uint64(initial_window_ms, 60000,
              "Initial window to guess a trace's rate");
DEFINE_double(decay_factor, 0.0, "How quickly to decay prediction");
DEFINE_double(link_capacity_scale, 1.0, "By how much to scale all links");
DEFINE_double(tm_scale, 1.0, "By how much to scale the traffic matrix");
DEFINE_string(opt, "CTR", "The optimizer to use");
DEFINE_double(generated_tm_scale, 1.2,
              "By how much to scale the generated matrix");
DEFINE_uint64(generated_tm_seed, 1ul, "Seed for the generated TM.");
DEFINE_double(generated_tm_max_global_utilization,
              std::numeric_limits<double>::max(),
              "Max global utilization for the generated matrix.");

static std::unique_ptr<nc::lp::DemandMatrix> GetMatrix(
    nc::net::GraphStorage* graph) {
  nc::lp::DemandGenerator generator(FLAGS_generated_tm_seed, graph);

  generator.AddUtilizationConstraint(0.1283987915407855, 0.0);
  generator.AddUtilizationConstraint(0.1521299093655589, 0.011384075954554769);
  generator.AddUtilizationConstraint(0.16832326283987917, 0.032535689078117526);
  generator.AddUtilizationConstraint(0.18342900302114803, 0.04879214954122174);
  generator.AddUtilizationConstraint(0.2006948640483384, 0.07481614717333394);
  generator.AddUtilizationConstraint(0.21255287009063445, 0.08781676191343549);
  generator.AddUtilizationConstraint(0.22658610271903323, 0.09920083786799026);
  generator.AddUtilizationConstraint(0.24169184290030213, 0.11058491382254503);
  generator.AddUtilizationConstraint(0.25463746223564954, 0.1219689897770998);
  generator.AddUtilizationConstraint(0.2632779456193353, 0.12848068122310513);
  generator.AddUtilizationConstraint(0.2805287009063444, 0.1398647571776599);
  generator.AddUtilizationConstraint(0.2934894259818731, 0.14960952619475876);
  generator.AddUtilizationConstraint(0.30534743202416914, 0.1561212176407641);
  generator.AddUtilizationConstraint(0.32045317220543806, 0.1642494478723162);
  generator.AddUtilizationConstraint(0.33555891238670693, 0.17563352382687097);
  generator.AddUtilizationConstraint(0.3474320241691843, 0.185401060995879);
  generator.AddUtilizationConstraint(0.3636102719033233, 0.19352929122743107);
  generator.AddUtilizationConstraint(0.37441087613293056, 0.19840167573598053);
  generator.AddUtilizationConstraint(0.38741691842900305, 0.20652990596753262);
  generator.AddUtilizationConstraint(0.39706948640483386, 0.21304159741353793);
  generator.AddUtilizationConstraint(0.4100151057401813, 0.22280913458254595);
  generator.AddUtilizationConstraint(0.41864048338368576, 0.2341932105371007);
  generator.AddUtilizationConstraint(0.4326737160120846, 0.24232144076865283);
  generator.AddUtilizationConstraint(0.44669184290030206, 0.25532205550875436);
  generator.AddUtilizationConstraint(0.4574924471299094, 0.2650895926777624);
  generator.AddUtilizationConstraint(0.4639577039274924, 0.268345438400765);
  generator.AddUtilizationConstraint(0.47690332326283985, 0.2829853600783225);
  generator.AddUtilizationConstraint(0.4909365558912387, 0.29111359030987455);
  generator.AddUtilizationConstraint(0.5017220543806646, 0.2976025136039708);
  generator.AddUtilizationConstraint(0.5146676737160121, 0.3089865895585256);
  generator.AddUtilizationConstraint(0.5254682779456193, 0.31549828100453087);
  generator.AddUtilizationConstraint(0.5330211480362538, 0.32688235695908563);
  generator.AddUtilizationConstraint(0.5448791540785498, 0.33664989412809365);
  generator.AddUtilizationConstraint(0.5535196374622356, 0.3447781243596458);
  generator.AddUtilizationConstraint(0.5621450151057401, 0.35290635459119785);
  generator.AddUtilizationConstraint(0.5718580060422961, 0.3610345848227499);
  generator.AddUtilizationConstraint(0.5837311178247734, 0.37567450650030737);
  generator.AddUtilizationConstraint(0.592356495468278, 0.38705858245486213);
  generator.AddUtilizationConstraint(0.6031570996978852, 0.39844265840941695);
  generator.AddUtilizationConstraint(0.6117824773413897, 0.4082101955784249);
  generator.AddUtilizationConstraint(0.6247280966767371, 0.42282734910407327);
  generator.AddUtilizationConstraint(0.6355135951661631, 0.4407231165046333);
  generator.AddUtilizationConstraint(0.649546827794562, 0.45210719245918807);
  generator.AddUtilizationConstraint(0.6614199395770393, 0.4683636529222923);
  generator.AddUtilizationConstraint(0.6700453172205438, 0.4781311900913003);
  generator.AddUtilizationConstraint(0.6840785498489427, 0.4943876505544045);
  generator.AddUtilizationConstraint(0.6948640483383686, 0.5074110334464152);
  generator.AddUtilizationConstraint(0.7088972809667674, 0.5285398784180688);
  generator.AddUtilizationConstraint(0.7218429003021148, 0.5431798000956263);
  generator.AddUtilizationConstraint(0.7326283987915407, 0.5626921062817332);
  generator.AddUtilizationConstraint(0.7531268882175226, 0.5887161039138453);
  generator.AddUtilizationConstraint(0.7682326283987916, 0.6082284100999522);
  generator.AddUtilizationConstraint(0.7768731117824773, 0.6277407162860591);
  generator.AddUtilizationConstraint(0.786570996978852, 0.6456364836866192);
  generator.AddUtilizationConstraint(0.8027643504531722, 0.6651487898727261);
  generator.AddUtilizationConstraint(0.8146374622356495, 0.6911727875048382);
  generator.AddUtilizationConstraint(0.8286555891238672, 0.7155802463514037);
  generator.AddUtilizationConstraint(0.8426888217522659, 0.7367090913230573);
  generator.AddUtilizationConstraint(0.8513141993957705, 0.7546048587236175);
  generator.AddUtilizationConstraint(0.8642598187311178, 0.7741171649097243);
  generator.AddUtilizationConstraint(0.8739728096676738, 0.7952687780332871);
  generator.AddUtilizationConstraint(0.885845921450151, 0.8098859315589354);
  generator.AddUtilizationConstraint(0.8944712990936555, 0.8277816989594955);
  generator.AddUtilizationConstraint(0.9031117824773414, 0.8489333120830583);
  generator.AddUtilizationConstraint(0.9138972809667674, 0.8700621570547119);
  generator.AddUtilizationConstraint(0.9225226586102719, 0.8895744632408188);
  generator.AddUtilizationConstraint(0.9322356495468278, 0.9058536918558321);
  generator.AddUtilizationConstraint(0.9387160120845921, 0.925365998041939);
  generator.AddUtilizationConstraint(0.9451812688821752, 0.9464948430135925);
  generator.AddUtilizationConstraint(0.9559818731117825, 0.962774071628606);
  generator.AddUtilizationConstraint(0.9635347432024169, 0.9790305320917102);
  generator.AddUtilizationConstraint(1.0, 0.9790305320917102);

  generator.AddHopCountLocalityConstraint(0.3, 2);
  generator.AddHopCountLocalityConstraint(0.2, 3);
  generator.AddHopCountLocalityConstraint(0.1, 4);
  generator.AddOutgoingFractionConstraint(1.0, 0.1 / FLAGS_generated_tm_scale);
  if (FLAGS_generated_tm_max_global_utilization !=
      std::numeric_limits<double>::max()) {
    generator.SetMaxGlobalUtilization(
        FLAGS_generated_tm_max_global_utilization);
  }

  generator.SetMinScaleFactor(1.3);
  std::unique_ptr<nc::lp::DemandMatrix> matrix =
      generator.GenerateMatrix(true, FLAGS_generated_tm_scale);
  CHECK(matrix) << "Unable to find matrix";
  return matrix;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  nc::metrics::InitMetrics();

  CHECK(!FLAGS_topology.empty());

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology), &node_order);
  builder.RemoveMultipleLinks();
  builder.ScaleCapacity(FLAGS_link_capacity_scale);
  nc::net::GraphStorage graph(builder);

  ctr::PcapTraceStore trace_store(FLAGS_pcap_trace_store);
  std::vector<ctr::BinSequence> all_bin_sequences;
  for (ctr::PcapDataTrace* trace : trace_store.AllTraces()) {
    all_bin_sequences.emplace_back(trace->ToSequence(trace->AllSlices()));
  }
  ctr::BinSequenceGenerator bin_sequence_generator(all_bin_sequences, 1000);

  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix;
  if (!FLAGS_traffic_matrix.empty()) {
    demand_matrix = nc::lp::DemandMatrix::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(FLAGS_traffic_matrix), node_order,
        &graph);
    demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);
  } else {
    demand_matrix = GetMatrix(&graph);
  }

  std::map<ctr::AggregateId, ctr::BinSequence> initial_sequences;
  for (const auto& matrix_element : demand_matrix->elements()) {
    ctr::BinSequence bin_sequence =
        ctr::BinsAtRate(matrix_element.demand,
                        std::chrono::milliseconds(FLAGS_initial_window_ms),
                        &bin_sequence_generator);
    LOG(ERROR) << bin_sequence.bin_count();

    initial_sequences.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(matrix_element.src, matrix_element.dst),
        std::forward_as_tuple(bin_sequence));
  }

  ctr::PathProvider path_provider(&graph);
  std::unique_ptr<ctr::Optimizer> opt;
  if (FLAGS_opt == "CTR") {
    opt = nc::make_unique<ctr::CTROptimizer>(&path_provider, false);
  } else if (FLAGS_opt == "B4") {
    opt = nc::make_unique<ctr::B4Optimizer>(&path_provider, false);
  } else if (FLAGS_opt == "B4(P)") {
    opt = nc::make_unique<ctr::B4Optimizer>(&path_provider, true);
  } else if (FLAGS_opt == "MinMax") {
    opt = nc::make_unique<ctr::MinMaxOptimizer>(&path_provider);
  }

  ctr::MeanScaleEstimatorFactory estimator_factory(
      {1.1, FLAGS_decay_factor, FLAGS_decay_factor, 10});
  ctr::RoutingSystem routing_system({}, opt.get(), &estimator_factory);

  ctr::NetMock net_mock(
      initial_sequences, std::chrono::milliseconds(FLAGS_period_duration_ms),
      std::chrono::milliseconds(FLAGS_history_bin_size_ms), &routing_system);
  net_mock.Run();
}
