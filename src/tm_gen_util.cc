#include <gflags/gflags.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/lp/mc_flow.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "common.h"
#include "opt/ctr.h"
#include "opt/opt.h"
#include "opt/path_provider.h"

DEFINE_string(topology, "", "A file with a topology");
DEFINE_string(output_suffix, "tm_gen",
              "Traffic matrices will be saved to files with this suffix");
DEFINE_double(scale, 1.2,
              "All aggregates in the generated matrix will initially fit on "
              "their shortest paths. The matrix can be further scaled after "
              "generation by this factor.");
DEFINE_double(min_scale_factor, 1.3,
              "The generated matrix should be scaleable by at least this much "
              "without becoming unfeasible.");
DEFINE_uint64(seed, 1ul, "Seed for the generated TM.");
DEFINE_uint64(max_try_count, 1000ul, "How many different TMs to try.");
DEFINE_uint64(tm_count, 1, "How many TMs to generate");
DEFINE_double(max_global_utilization, std::numeric_limits<double>::max(),
              "Max global utilization for the generated matrix.");
DEFINE_double(link_capacity_scale, 1.0,
              "By how much to scale all links' capacity");
DEFINE_uint64(max_variable_count, 10000, "Max number of LP variables");
DEFINE_uint64(max_constraint_count, 10000, "Max number of LP constraints");

// Returns the cost of a matrix. The cost is tied to the improvement expected
// from using CTR vs MinMax.
static double Cost(const nc::lp::DemandMatrix& matrix) {
  nc::net::Bandwidth max_element = nc::net::Bandwidth::Zero();
  for (const auto& element : matrix.elements()) {
    max_element = std::max(element.demand, max_element);
  }

  std::map<nc::lp::SrcAndDst, size_t> flow_counts;
  for (const auto& element : matrix.elements()) {
    double f = element.demand.Mbps() / max_element.Mbps();
    size_t flow_count = f * 1000;
    flow_counts[{element.src, element.dst}] = std::max(1ul, flow_count);
  }

  ctr::TrafficMatrix tm(matrix, flow_counts);
  ctr::PathProvider path_provider(tm.graph());
  ctr::CTROptimizer ctr_optimizer(&path_provider, 1.0, false, false);
  std::unique_ptr<ctr::RoutingConfiguration> ctr_out =
      ctr_optimizer.Optimize(tm);

  ctr::MinMaxOptimizer minmax_optimizer(&path_provider, 1.0, false);
  std::unique_ptr<ctr::RoutingConfiguration> minmax_out =
      minmax_optimizer.Optimize(tm);

  ctr::RoutingConfigurationDelta delta = ctr_out->GetDifference(*minmax_out);
  return delta.total_per_flow_delay_delta;
}

static std::vector<std::unique_ptr<nc::lp::DemandMatrix>> GetMatrices(
    nc::net::GraphStorage* graph) {
  nc::lp::DemandGenerator generator(FLAGS_seed, graph);

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
  generator.AddOutgoingFractionConstraint(1.0, 0.1 / FLAGS_scale);
  if (FLAGS_max_global_utilization != std::numeric_limits<double>::max()) {
    generator.SetMaxGlobalUtilization(FLAGS_max_global_utilization);
  }

  generator.SetMinScaleFactor(FLAGS_min_scale_factor);
  generator.SetMinOverloadedLinkCount(2);

  size_t var_count, constraint_count;
  std::tie(var_count, constraint_count) = generator.EstimateLPSize();
  LOG(INFO) << "Estimated LP size " << var_count << " / " << constraint_count;
  if (var_count > FLAGS_max_variable_count ||
      constraint_count > FLAGS_max_constraint_count) {
    LOG(FATAL) << "LP too large";
  }

  std::vector<std::unique_ptr<nc::lp::DemandMatrix>> matrices =
      generator.GenerateMatrices(FLAGS_max_try_count, FLAGS_scale);

  std::sort(matrices.begin(), matrices.end(),
            [](const std::unique_ptr<nc::lp::DemandMatrix>& lhs,
               const std::unique_ptr<nc::lp::DemandMatrix>& rhs) {
              return lhs->SPGlobalUtilization() > rhs->SPGlobalUtilization();
            });

  matrices.resize(
      std::min(static_cast<size_t>(FLAGS_tm_count), matrices.size()));
  return matrices;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_topology.empty());

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology), &node_order);
  builder.RemoveMultipleLinks();
  builder.ScaleCapacity(FLAGS_link_capacity_scale);
  nc::net::GraphStorage graph(builder);
  std::vector<std::unique_ptr<nc::lp::DemandMatrix>> demand_matrices =
      GetMatrices(&graph);
  for (uint64_t i = 0; i < demand_matrices.size(); ++i) {
    const nc::lp::DemandMatrix& demand_matrix = *demand_matrices[i];
    std::string output = nc::StrCat(FLAGS_topology, "_", FLAGS_output_suffix,
                                    "_", i, ".demands");

    nc::File::WriteStringToFileOrDie(demand_matrix.ToRepetita(node_order),
                                     output);
    LOG(INFO) << "Written TM to " << output;
  }
}
