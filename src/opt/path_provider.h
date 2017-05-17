#ifndef CTR_PATH_PROVIDER_H
#define CTR_PATH_PROVIDER_H

#include "ncode_net/src/net_common.h"
#include "ncode_net/src/graph_query.h"

namespace ctr {

// Can provide the K next shortest paths on demand for an aggregate.
class PathProvider {
 public:
  // Returns the K shortest paths computed so far for an aggregate. Empty if
  // NextShortestPathOrNull has not been called yet.
  const std::vector<const nc::net::Walk*>& KShorestPaths(AggregateId aggregate);

  // Extends the K shortest paths for an aggregate and returns the path.
  const nc::net::Walk* NextShortestPathOrNull(AggregateId aggregate);

  // Returns a path that avoids the given set of links.
  const nc::net::Walk* AvoidingPathOrNull(
      AggregateId aggregate, const nc::net::GraphLinkSet& to_avoid);

  // Takes ownership of a path.
  const nc::net::Walk* TakeOwnership(std::unique_ptr<nc::net::Walk> path);
};

}  // namespace ctr

#endif
