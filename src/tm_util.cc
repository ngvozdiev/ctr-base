#include <gflags/gflags.h>
#include <memory>
#include <string>
#include <vector>

#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"

DEFINE_string(topology, "", "A file with a topology");
DEFINE_string(traffic_matrix, "", "A file with a traffic matrix");
DEFINE_double(max_demand_to_max_link_capacity, 0.5,
              "The fraction max_demand / max_link_capacity");
DEFINE_double(
    target_commodity_scale_factor, 1.0,
    "The matrix should be feasible if all commodities are scaled by this much");
DEFINE_uint64(seed, 1, "Seed for the RNG");
DEFINE_string(output, "", "If non-empty will save the output matrix there");

static void MakeFitRecursive(
    const nc::lp::DemandMatrix& input,
    const std::vector<nc::lp::DemandMatrixElement>& all_elements, size_t from,
    size_t to, double target_csf,
    std::unique_ptr<nc::lp::DemandMatrix>* best_so_far, double* best_csf) {
  size_t range = (to - from) / 2;
  size_t threshold = from + range;
  if (threshold == 0) {
    return;
  }

  if (threshold == all_elements.size()) {
    return;
  }

  std::set<nc::lp::DemandMatrix::NodePair> to_remove;
  for (size_t i = 0; i < threshold; ++i) {
    const nc::lp::DemandMatrixElement& element = all_elements[i];
    to_remove.emplace(element.src, element.dst);
  }

  std::unique_ptr<nc::lp::DemandMatrix> candidate =
      input.RemovePairs(to_remove);
  if (range == 0) {
    return;
  }

  if (candidate->IsFeasible({})) {
    double csf = candidate->MaxCommodityScaleFractor();
    if (std::abs(*best_csf - target_csf) > std::abs(csf - target_csf)) {
      *best_so_far = std::move(candidate);
      *best_csf = csf;
    }

    if (std::abs(csf - target_csf) < 0.01) {
      return;
    }

    if (csf > target_csf) {
      MakeFitRecursive(input, all_elements, from, threshold, target_csf,
                       best_so_far, best_csf);
    } else {
      MakeFitRecursive(input, all_elements, threshold, to, target_csf,
                       best_so_far, best_csf);
    }

    return;
  }

  MakeFitRecursive(input, all_elements, threshold, to, target_csf, best_so_far,
                   best_csf);
}

// Removes small commodities until the TM fits with a given
// commodity_scale_factor.
static std::unique_ptr<nc::lp::DemandMatrix> MakeFit(
    const nc::lp::DemandMatrix& input, double target_csf) {
  // Demands, ordered by size.
  std::vector<nc::lp::DemandMatrixElement> all_elements = input.elements();

  std::mt19937 rnd(FLAGS_seed);
  std::shuffle(all_elements.begin(),
               std::next(all_elements.begin(), all_elements.size()), rnd);

  nc::net::Bandwidth max_demand;
  for (const auto& element : all_elements) {
    max_demand = std::max(max_demand, element.demand);
  }

  nc::net::Bandwidth max_link_capacity;
  const nc::net::GraphStorage* graph = input.graph();
  for (nc::net::GraphLinkIndex link_index : graph->AllLinks()) {
    max_link_capacity =
        std::max(max_link_capacity, graph->GetLink(link_index)->bandwidth());
  }

  double scale =
      FLAGS_max_demand_to_max_link_capacity * (max_link_capacity / max_demand);
  auto scaled_input = input.Scale(scale);

  std::unique_ptr<nc::lp::DemandMatrix> best_candidate;
  double best_csf = std::numeric_limits<double>::max();
  MakeFitRecursive(*scaled_input, all_elements, 0, all_elements.size(),
                   target_csf, &best_candidate, &best_csf);
  return best_candidate;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::vector<std::string> node_order;
  nc::net::GraphBuilder builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology), &node_order);
  builder.RemoveMultipleLinks();

  nc::net::GraphStorage graph(builder);
  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix =
      nc::lp::DemandMatrix::LoadRepetitaOrDie(
          nc::File::ReadFileToStringOrDie(FLAGS_traffic_matrix), node_order,
          &graph);
  demand_matrix = MakeFit(*demand_matrix, FLAGS_target_commodity_scale_factor);

  LOG(INFO) << graph.Stats().ToString();
  if (demand_matrix) {
    LOG(INFO) << demand_matrix->ToString();
  }

  if (!FLAGS_output.empty()) {
    std::string serialized_matrix = demand_matrix->ToRepetita(node_order);
    nc::File::WriteStringToFile(serialized_matrix, FLAGS_output);
  }
}
