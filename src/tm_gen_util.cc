#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/algorithm.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/strutil.h"
#include "ncode_common/src/viz/grapher.h"

DEFINE_string(topology, "", "A file with a topology");
DEFINE_string(output_suffix, "tm_gen",
              "Traffic matrices will be saved to files with this suffix");
DEFINE_double(min_scale_factor, 1.3,
              "The generated matrix should be scaleable by at least this much "
              "without becoming unfeasible.");
DEFINE_double(locality, 0.0,
              "How much of does locality factor into the traffic matrix.");
DEFINE_uint64(seed, 1ul, "Seed for the generated TM.");
DEFINE_uint64(tm_count, 1, "How many TMs to generate");

static void PlotRateVsDistance(const nc::lp::DemandMatrix& matrix) {
  const nc::net::GraphStorage* graph = matrix.graph();
  nc::net::AllPairShortestPath sp({}, graph->AdjacencyList(), nullptr, nullptr);

  std::vector<nc::lp::DemandMatrixElement> elements = matrix.elements();
  std::sort(elements.begin(), elements.end(),
            [&sp](const nc::lp::DemandMatrixElement& lhs,
                  const nc::lp::DemandMatrixElement& rhs) {
              nc::net::Delay lhs_distance = sp.GetDistance(lhs.src, lhs.dst);
              nc::net::Delay rhs_distance = sp.GetDistance(rhs.src, rhs.dst);
              return lhs_distance < rhs_distance;
            });

  nc::viz::DataSeries2D to_plot;
  double cumulative_rate = 0;
  for (const auto& element : elements) {
    nc::net::Delay distance = sp.GetDistance(element.src, element.dst);
    cumulative_rate += element.demand.Mbps();
    to_plot.data.emplace_back(distance.count(), cumulative_rate);
  }

  nc::viz::PythonGrapher grapher("rate_vs_distance");
  grapher.PlotLine({}, {to_plot});
}

static void PlotAggregateRates(const nc::lp::DemandMatrix& matrix) {
  std::vector<double> to_plot;
  for (const auto& element : matrix.elements()) {
    to_plot.emplace_back(element.demand.Mbps());
  }

  nc::viz::PythonGrapher grapher("aggregate_rates");
  grapher.PlotCDF({}, {{"", to_plot}});
}

static std::vector<std::unique_ptr<nc::lp::DemandMatrix>> GetMatrices(
    nc::net::GraphStorage* graph) {
  nc::lp::DemandGenerator generator(graph, FLAGS_seed);

  std::vector<std::unique_ptr<nc::lp::DemandMatrix>> matrices;
  for (size_t i = 0; i < FLAGS_tm_count; ++i) {
    auto tm = generator.Generate(FLAGS_min_scale_factor, FLAGS_locality);

    // Will ignore all aggregates less than 1Mbps.
    tm = tm->Filter([](const nc::lp::DemandMatrixElement& element) {
      return element.demand < nc::net::Bandwidth::FromMBitsPerSecond(1);
    });

    matrices.emplace_back(std::move(tm));
  }

  return matrices;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_topology.empty());

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology), &node_order);
  builder.RemoveMultipleLinks();
  nc::net::GraphStorage graph(builder);
  std::vector<std::unique_ptr<nc::lp::DemandMatrix>> demand_matrices =
      GetMatrices(&graph);
  for (uint64_t i = 0; i < demand_matrices.size(); ++i) {
    const nc::lp::DemandMatrix& demand_matrix = *demand_matrices[i];
    std::string output = nc::StrCat(FLAGS_topology, "_", FLAGS_output_suffix,
                                    "_", i, ".demands");

    nc::File::WriteStringToFileOrDie(demand_matrix.ToRepetita(node_order),
                                     output);
    PlotRateVsDistance(demand_matrix);
    PlotAggregateRates(demand_matrix);

    LOG(INFO) << "TM " << demand_matrix.ToString();
    LOG(INFO) << "Written TM to " << output;
  }
}
