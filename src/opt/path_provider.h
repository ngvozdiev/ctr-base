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

  // Adds paths in K shortest order until a path is added that avoids a set of
  // links or max_count is reached. Paths are added starting at start_k.
  std::vector<const nc::net::Walk*> KShortestUntilAvoidingPath(
      const AggregateId& aggregate, const nc::net::GraphLinkSet& to_avoid,
      size_t start_k, size_t max_count) {
    std::vector<const nc::net::Walk*> out;

    Generator* generator = FindOrCreateGenerator(aggregate);
    size_t i = start_k;
    while (out.size() != max_count) {
      const nc::net::Walk* next_path = generator->KthShortestPathOrNull(i++);
      if (next_path == nullptr) {
        return {};
      }

      out.emplace_back(next_path);
      if (next_path->ContainsAny(to_avoid)) {
        continue;
      }

      break;
    }

    return out;
  }

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
