#ifndef CTR_OPT_COMPARE_H
#define CTR_OPT_COMPARE_H

#include "../common.h"
#include "opt.h"

namespace ctr {

// Returns a single line with the following space-separated information:
// number of aggregates
// distribution of absolute path stretches in ms (100 percentile values)
// distribution of relative path stretches
// distribution of numbers of paths per aggregate
// distribution of link utilization
// distribution of unmet demand per link
std::string DumpRoutingInfo(const RoutingConfiguration& routing);

// Returns a single line with space-separated information comparing two routing
// configurations.
// fraction of total volume that changed paths
// fraction of total flow count that changed paths
// distribution of per-aggregate fraction changes
// number of route adds
// number of route updates
// number of route removals
std::string DumpRoutingDelta(const RoutingConfiguration& routing_one,
                             const RoutingConfiguration& routing_two);

}  // namespace ctr

#endif
