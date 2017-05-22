#include "path_provider.h"

#include <algorithm>
#include <set>

#include "ncode_common/src/common.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/net/algorithm.h"

namespace ctr {

PathProvider::Generator* PathProvider::FindOrCreateGenerator(
    const AggregateId& aggregate_id) {
  PathProvider::Generator* generator_ptr =
      nc::FindSmartPtrOrNull(generators_, aggregate_id);
  if (generator_ptr != nullptr) {
    return generator_ptr;
  }

  std::set<nc::net::ConstraintSet> constraint_sets;
  nc::net::ConstraintSet empty_set;
  constraint_sets.emplace(empty_set);

  auto new_gen = nc::make_unique<PathProvider::Generator>(
      aggregate_id.src(), aggregate_id.dst(), *graph_, constraint_sets);
  PathProvider::Generator* raw_ptr = new_gen.get();

  generators_[aggregate_id] = std::move(new_gen);
  return raw_ptr;
}

const std::vector<const nc::net::Walk*>& PathProvider::KShorestPaths(
    const AggregateId& aggregate) {
  PathProvider::Generator* generator = FindOrCreateGenerator(aggregate);
  return generator->k_paths();
}

const nc::net::Walk* PathProvider::AvoidingPathOrNull(
    const AggregateId& aggregate, const nc::net::GraphLinkSet& to_avoid) {
  PathProvider::Generator* generator = FindOrCreateGenerator(aggregate);
  auto path = generator->ShortestPathThatAvoids({}, to_avoid);
  if (!path) {
    return nullptr;
  }

  return TakeOwnership(std::move(path));
}

const nc::net::Walk* PathProvider::TakeOwnership(
    std::unique_ptr<nc::net::Walk> path) {
  auto it = owned_paths_.find(path);
  if (it != owned_paths_.end()) {
    return it->get();
  }

  const nc::net::Walk* raw_ptr = path.get();
  owned_paths_.insert(std::move(path));
  return raw_ptr;
}

}  // namespace ctr
