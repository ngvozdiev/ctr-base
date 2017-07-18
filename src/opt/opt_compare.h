#ifndef CTR_OPT_COMPARE_H
#define CTR_OPT_COMPARE_H

#include "ncode_common/src/viz/grapher.h"

namespace ctr {
class RoutingConfiguration;
struct RoutingConfigurationDelta;
} /* namespace ctr */

namespace ctr {

// Records various statistics about a RoutingConfiguration to metrics.
void RecordRoutingConfig(const std::string& topology, const std::string& tm,
                         const std::string& opt,
                         const RoutingConfiguration& routing);

// Information about the difference between two routing configs. Has fields:
// fraction of total volume that changed paths
// fraction of total flow count that changed paths
// number of route adds
// number of route updates
// number of route removals
// distribution of per-aggregate fraction changes
class RoutingConfigDeltaInfo : public nc::viz::NpyArray {
 public:
  RoutingConfigDeltaInfo();

  void Add(const RoutingConfigurationDelta& delta);
};

}  // namespace ctr

#endif
