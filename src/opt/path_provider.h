#ifndef CTR_PATH_PROVIDER_H
#define CTR_PATH_PROVIDER_H

#include "ncode_net/src/net_common.h"
#include "ncode_net/src/graph_query.h"

namespace ctr {

// Can provide the K next shortest paths on demand for an aggregate.
class PathProvider {
 public:
  // Returns the K shortest paths. The returned vector will have up to k + 1
  // elements. K=0 is the shortest path.
  std::vector<const nc::net::Walk*> KShorestPaths(AggregateId aggregate,
                                                  size_t k);
};

}  // namespace ctr

#endif
