#ifndef CTR_PATH_PROVIDER_H
#define CTR_PATH_PROVIDER_H

#include <map>
#include <memory>
#include <vector>

#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/algorithm.h"
#include "../common.h"

namespace ctr {

// Can provide the K next shortest paths on demand for an aggregate.
class PathProvider {
 public:
  PathProvider(const nc::net::GraphStorage* graph) : graph_(graph) {}

  // Returns the K shortest paths computed so far for an aggregate. Empty if
  // NextShortestPathOrNull has not been called yet.
  const std::vector<const nc::net::Walk*>& KShorestPaths(
      const AggregateId& aggregate);

  // Extends the K shortest paths for an aggregate and returns the path.
  const nc::net::Walk* NextShortestPathOrNull(const AggregateId& aggregate);

  // Returns a path that avoids the given set of links.
  const nc::net::Walk* AvoidingPathOrNull(
      const AggregateId& aggregate, const nc::net::GraphLinkSet& to_avoid);

  // Takes ownership of a path.
  const nc::net::Walk* TakeOwnership(std::unique_ptr<nc::net::Walk> path);

  const nc::net::GraphStorage* graph() const { return graph_; };

 private:
  using Generator = nc::net::DisjunctKShortestPathsGenerator;

  Generator* FindOrCreateGenerator(const AggregateId& aggregate_id);

  // The graph.
  const nc::net::GraphStorage* graph_;

  // A path generator for each aggregate.
  std::map<AggregateId, std::unique_ptr<Generator>> generators_;

  // Paths owned by this object should be compared by their value, not by the
  // value of their pointers.
  struct OwnedPathsComparator {
    bool operator()(const std::unique_ptr<nc::net::Walk>& lhs,
                    const std::unique_ptr<nc::net::Walk>& rhs) const {
      return *lhs < *rhs;
    }
  };

  // A set of paths owned by this object.
  std::set<std::unique_ptr<nc::net::Walk>, OwnedPathsComparator> owned_paths_;

  DISALLOW_COPY_AND_ASSIGN(PathProvider);
};

}  // namespace ctr

#endif
