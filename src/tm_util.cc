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
DEFINE_double(tm_scale, 1.0, "The traffic matrix will be scaled by this much");
DEFINE_double(
    target_commodity_scale_factor, 1.0,
    "The matrix should be feasible if all commodities are scaled by this much");
DEFINE_uint64(seed, 1, "Seed for the RNG");

static std::unique_ptr<nc::lp::DemandMatrix> MakeFitRecursive(
    const nc::lp::DemandMatrix& input,
    const std::vector<nc::lp::DemandMatrixElement>& all_elements, size_t from,
    size_t to, double target_commodity_scale_factor) {
  size_t range = (to - from) / 2;
  size_t threshold = from + range;
  if (threshold == 0) {
    return {};
  }

  if (threshold == all_elements.size()) {
    return {};
  }

  std::set<nc::lp::DemandMatrix::NodePair> to_remove;
  for (size_t i = 0; i < threshold; ++i) {
    const nc::lp::DemandMatrixElement& element = all_elements[i];
    to_remove.emplace(element.src, element.dst);
  }

  std::unique_ptr<nc::lp::DemandMatrix> candidate =
      input.RemovePairs(to_remove);
  if (range == 0) {
    return candidate;
  }

  if (candidate->IsFeasible({})) {
    double sf = candidate->MaxCommodityScaleFractor();
    LOG(ERROR) << "T " << threshold << "/" << all_elements.size() << " sf "
               << sf;
    if (std::abs(sf - target_commodity_scale_factor) < 0.01) {
      return candidate;
    }

    if (sf > target_commodity_scale_factor) {
      return MakeFitRecursive(input, all_elements, from, threshold,
                              target_commodity_scale_factor);
    } else {
      return MakeFitRecursive(input, all_elements, threshold, to,
                              target_commodity_scale_factor);
    }
  }

  LOG(ERROR) << "T " << threshold << "/" << all_elements.size() << " NF";
  return MakeFitRecursive(input, all_elements, threshold, to,
                          target_commodity_scale_factor);
}

// Removes small commodities until the TM fits with a given
// commodity_scale_factor.
static std::unique_ptr<nc::lp::DemandMatrix> MakeFit(
    const nc::lp::DemandMatrix& input, double target_commodity_scale_factor) {
  // Demands, ordered by size.
  std::vector<nc::lp::DemandMatrixElement> all_elements = input.elements();

  std::mt19937 rnd(FLAGS_seed);
  std::shuffle(all_elements.begin(),
               std::next(all_elements.begin(), all_elements.size()), rnd);

  return MakeFitRecursive(input, all_elements, 0, all_elements.size(),
                          target_commodity_scale_factor);
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
  demand_matrix = demand_matrix->Scale(FLAGS_tm_scale);
  demand_matrix = MakeFit(*demand_matrix, FLAGS_target_commodity_scale_factor);

  LOG(INFO) << graph.Stats().ToString();
  if (demand_matrix) {
    LOG(INFO) << demand_matrix->ToString();
  }
}
