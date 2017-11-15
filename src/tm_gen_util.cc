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

DEFINE_string(topology_files, "", "Topology files");
DEFINE_string(output_pattern,
              "demand_matrices/scale_factor_$2/locality_$1/$0_$3.demands",
              "Traffic matrices will be saved to files named after this "
              "pattern, with $0 replaced by the topology name, $1 replaced by "
              "locality, $2 replaced by scale factor and $3 replaced by a "
              "unique integer identifier.");
DEFINE_double(min_scale_factor, 1.3,
              "The generated matrix should be scaleable by at least this much "
              "without becoming unfeasible.");
DEFINE_double(locality, 0.0,
              "How much of does locality factor into the traffic matrix.");
DEFINE_uint64(seed, 1ul, "Seed for the generated TM.");
DEFINE_uint64(tm_count, 1, "Number of TMs to generate.");
DEFINE_uint64(threads, 2, "Number of threads to use.");

using DemandMatrixVector = std::vector<std::unique_ptr<nc::lp::DemandMatrix>>;

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

static DemandMatrixVector GetMatrices(
    const nc::lp::DemandGenerator& demand_generator, size_t count,
    uint64_t seed) {
  std::mt19937 gen(seed);
  std::vector<std::unique_ptr<nc::lp::DemandMatrix>> matrices;
  for (size_t i = 0; i < count; ++i) {
    auto tm =
        demand_generator.Generate(FLAGS_min_scale_factor, FLAGS_locality, &gen);

    // Will ignore all aggregates less than 1Mbps.
    tm = tm->Filter([](const nc::lp::DemandMatrixElement& element) {
      return element.demand < nc::net::Bandwidth::FromMBitsPerSecond(1);
    });

    matrices.emplace_back(std::move(tm));
  }

  return matrices;
}

static DemandMatrixVector GetAllMatrices(const nc::net::GraphStorage* graph) {
  nc::lp::DemandGenerator demand_generator;
  std::vector<uint64_t> seeds;
  for (size_t i = 0; i < FLAGS_threads; ++i) {
    seeds.emplace_back(FLAGS_seed + i);
  }

  // This may cause the total to not be equal to FLAGS_tm_count.
  size_t per_thread_count = FLAGS_tm_count / FLAGS_threads;

  std::vector<std::unique_ptr<DemandMatrixVector>> all_matrices(FLAGS_threads);
  nc::RunInParallelWithResult(
      seeds, [&demand_generator, per_thread_count](uint64_t seed) {
        DemandMatrixVector matrices_for_thread =
            GetMatrices(demand_generator, per_thread_count, seed);
        return nc::make_unique<DemandMatrixVector>(matrices_for_thread);
      });

  // Time to merge the matrices into a final output vector.
  DemandMatrixVector out;
  for (auto& matrices_for_thread : all_matrices) {
    DemandMatrixVector& vector_for_thread = *matrices_for_thread;
    out.insert(out.end(), vector_for_thread.begin(), vector_for_thread.end());
  }

  return out;
}

// Globs filenames from FLAGS_topology_files
static std::vector<std::string> GetTopologyFiles() {
  std::vector<std::string> out;
  std::vector<std::string> split = nc::Split(FLAGS_topology_files, ",");
  for (const std::string& piece : split) {
    std::vector<std::string> files = nc::Glob(piece);
    out.insert(out.end(), files.begin(), files.end());
  }

  return out;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::vector<std::string> topology_files = GetTopologyFiles();
  CHECK(!topology_files.empty());

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology), &node_order);
  builder.RemoveMultipleLinks();
  nc::net::GraphStorage graph(builder);
  std::vector<std::unique_ptr<nc::lp::DemandMatrix>> demand_matrices =
      GetAllMatrices(&graph);
  for (uint64_t i = 0; i < demand_matrices.size(); ++i) {
    const nc::lp::DemandMatrix& demand_matrix = *demand_matrices[i];
    nc::Substitute(FLAGS_output_pattern)

        std::string output = nc::StrCat(
            FLAGS_topology, "_", FLAGS_output_suffix, "_", i, ".demands");

    nc::File::WriteStringToFileOrDie(demand_matrix.ToRepetita(node_order),
                                     output);
    PlotRateVsDistance(demand_matrix);
    PlotAggregateRates(demand_matrix);

    LOG(INFO) << "TM " << demand_matrix.ToString();
    LOG(INFO) << "Written TM to " << output;
  }
}
