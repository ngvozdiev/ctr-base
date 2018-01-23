#include <gflags/gflags.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/lp/demand_matrix.h>
#include <ncode/map_util.h>
#include <ncode/net/net_common.h>
#include <ncode/strutil.h>
#include <ncode/substitute.h>
#include <chrono>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "common.h"
#include "demand_matrix_input.h"
#include "opt/opt.h"
#include "opt/path_provider.h"
#include "topology_input.h"

DEFINE_double(sp_fraction, 1.2, "How far from the SP a path can be");
DEFINE_double(link_fraction_limit, 0.8,
              "At least this much of the SP's links can be routed around");
DEFINE_string(opt, "SP,B4,MinMaxK10,CTR,MinMaxLD", "The optimizers to plot");
DEFINE_string(output, "", "The file to store data to");

struct RCSummary {
  RCSummary(double change_in_delay, double fraction_congested,
            double routability, const std::string* topology_file,
            const std::string* demand_file)
      : change_in_delay(change_in_delay),
        fraction_congested(fraction_congested),
        topology_file(topology_file),
        demand_file(demand_file),
        topology_routability(routability) {}

  double change_in_delay;
  double fraction_congested;

  // Link utilizations.
  nc::net::GraphLinkMap<double> utilizations;

  // The topology/demand matrix.
  const std::string* topology_file;
  const std::string* demand_file;

  // The routability metric of the topology.
  double topology_routability;
};

struct Input {
  const ctr::DemandMatrixAndFilename* demand;
  const ctr::TopologyAndFilename* topology;
  double routability;
};

// Will return the change in total delay.
static double GetDelayDatapointForRC(const ctr::RoutingConfiguration& rc) {
  nc::net::Delay total = rc.TotalPerFlowDelay(false);
  nc::net::Delay total_sp = rc.TotalPerFlowDelay(true);
  double fraction = static_cast<double>(total.count()) / total_sp.count();
  return fraction;
}

// Will return the fraction of aggregates that do not fit.
static double GetCapacityDatapointForRC(const ctr::RoutingConfiguration& rc) {
  double total_count = rc.demands().size();
  return rc.OverloadedAggregates() / total_count;
}

static double GetDatapointForTopology(const nc::net::GraphStorage& graph) {
  return ctr::GetFractionOfPairsAboveLinkFraction(graph, FLAGS_sp_fraction,
                                                  FLAGS_link_fraction_limit);
}

static std::string GetFilename(const std::string& tm_file,
                               const std::string opt_string) {
  std::string tm_base = nc::StringReplace(tm_file, ".demands", "", true);
  return nc::StrCat(tm_base, "_", opt_string, ".rc");
}

static std::vector<RCSummary> ParseRcs(const std::string& opt,
                                       const std::vector<Input>& matrices) {
  std::vector<RCSummary> out;
  for (const auto& input : matrices) {
    const nc::lp::DemandMatrix* demand_matrix =
        input.demand->demand_matrix.get();
    const nc::net::GraphStorage* graph = demand_matrix->graph();
    const std::string& topology_file = input.topology->file;
    const std::string& demand_file = input.demand->file;
    const std::vector<std::string>& node_order = input.topology->node_order;

    ctr::PathProvider path_provider(graph);
    std::string rc_filename = GetFilename(input.demand->file, opt);
    if (!nc::File::Exists(rc_filename)) {
      LOG(INFO) << "Missing " << rc_filename;
      continue;
    }

    std::string rc_serialized = nc::File::ReadFileToStringOrDie(rc_filename);
    std::unique_ptr<ctr::TrafficMatrix> tm =
        ctr::TrafficMatrix::DistributeFromDemandMatrix(*demand_matrix);
    LOG(INFO) << "Will parse " << demand_file << " at " << topology_file
              << " opt " << opt;
    auto rc = ctr::RoutingConfiguration::LoadFromSerializedText(
        *tm, node_order, rc_serialized, &path_provider);

    double delay_datapoint = GetDelayDatapointForRC(*rc);
    double capacity_datapoint = GetCapacityDatapointForRC(*rc);

    out.emplace_back(delay_datapoint, capacity_datapoint, input.routability,
                     &topology_file, &demand_file);
  }

  return out;
}

static std::string Dump(const std::string& opt,
                        const std::vector<RCSummary>& rcs) {
  std::string out = "";
  for (const auto& rc : rcs) {
    std::string topology_file = *(rc.topology_file);
    std::string demand_file = *(rc.demand_file);
    nc::SubstituteAndAppend(&out, "$0 $1 $2 $3 $4 $5\n", opt,
                            rc.topology_routability, rc.change_in_delay,
                            rc.fraction_congested, topology_file, demand_file);
  }

  return out;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(FLAGS_output != "") << "Need output";

  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> matrices;
  std::tie(topologies, matrices) = ctr::GetDemandMatrixInputs(false);

  std::map<std::string, const ctr::TopologyAndFilename*> topologies_by_name;
  std::map<std::string, double> topology_datapoints;
  for (const auto& topology : topologies) {
    topologies_by_name[topology.file] = &topology;

    LOG(INFO) << "Getting datapoint for " << topology.file;
    double datapoint = GetDatapointForTopology(*topology.graph);
    topology_datapoints[topology.file] = datapoint;
  }

  std::vector<Input> inputs;
  for (const ctr::DemandMatrixAndFilename& matrix : matrices) {
    const ctr::TopologyAndFilename* topology =
        nc::FindOrDie(topologies_by_name, matrix.topology_file);
    double routability =
        nc::FindOrDie(topology_datapoints, matrix.topology_file);
    inputs.push_back({&matrix, topology, routability});
  }

  std::string total = "";
  std::vector<std::string> opts = nc::Split(FLAGS_opt, ",", true);
  for (const std::string& opt : opts) {
    std::vector<RCSummary> rcs = ParseRcs(opt, inputs);
    total += Dump(opt, rcs);
  }

  nc::File::WriteStringToFileOrDie(total, FLAGS_output);
}
