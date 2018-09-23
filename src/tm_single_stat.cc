#include <gflags/gflags.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ncode/common.h"
#include "ncode/file.h"
#include "ncode/lp/demand_matrix.h"
#include "ncode/net/algorithm.h"
#include "ncode/net/net_common.h"
#include "ncode/thread_runner.h"
#include "ncode/viz/ctemplate/template.h"
#include "ncode/viz/ctemplate/template_dictionary.h"
#include "ncode/viz/ctemplate/template_enums.h"
#include "ncode/viz/grapher.h"
#include "demand_matrix_input.h"
#include "topology_input.h"
#include "plot_algorithm_eval_tools.h"

DEFINE_string(opt_eval_metrics, "", "Location of metrics from opt_eval_util");
DEFINE_string(optimizers, "MinMaxLD,CTRNFC,B4,CTR,MinMaxK10",
              "Optimizers to plot.");

using namespace ctr::alg_eval;

static void PlotSingleTMStats(const ctr::DemandMatrixAndFilename& input,
                              const DataStorage& data_storage,
                              const std::string& id) {
  const std::map<std::string, TMStateMap>& data = data_storage.data();
  nc::viz::LinePlot path_stretch_plot =
      GetLinePlot("milliseconds", "CDF", "Path stretch (absolute)");

  std::vector<std::string> optimizers = nc::Split(FLAGS_optimizers, ",");
  for (const std::string& opt : optimizers) {
    const TMStateMap& state_map = nc::FindOrDie(data, opt);
    TopologyAndTM key(input.topology_file, input.file);

    const OptimizerTMState& tm_state = *(nc::FindOrDieNoPrint(state_map, key));
    std::vector<double> percentiles = GetStretchDistribution(tm_state);
    std::vector<std::pair<double, double>> values;
    for (size_t i = 0; i < percentiles.size(); ++i) {
      double p = percentiles[i];
      values.emplace_back(p, i);
    }
    path_stretch_plot.AddData(opt, values);
  }
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::vector<ctr::TopologyAndFilename> topologies;
  std::vector<ctr::DemandMatrixAndFilename> demands;

  CHECK(demands.size() == 1);
  PlotDemandStats(demands);
}
