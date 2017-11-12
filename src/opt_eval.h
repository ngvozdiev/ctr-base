#ifndef CTR_OPT_EVAL_H
#define CTR_OPT_EVAL_H

#include <iostream>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_common.h"

namespace ctr {

struct OptEvalInput {
  OptEvalInput(const std::string& topology_file, const std::string& tm_file,
               std::unique_ptr<nc::lp::DemandMatrix> demand_matrix)
      : topology_file(topology_file),
        tm_file(tm_file),
        demand_matrix(std::move(demand_matrix)) {}

  OptEvalInput() {}

  std::string topology_file;
  std::string tm_file;
  std::unique_ptr<nc::lp::DemandMatrix> demand_matrix;
};

std::pair<std::vector<std::unique_ptr<nc::net::GraphStorage>>,
          std::vector<OptEvalInput>>
GetOptEvalInputs();

}  // namespace ctr

#endif
