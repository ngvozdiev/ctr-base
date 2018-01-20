#ifndef CTR_COMMON_H
#define CTR_COMMON_H

#include <ncode/common.h>
#include <ncode/htsim/match.h>
#include <ncode/htsim/packet.h>
#include <ncode/lp/mc_flow.h>
#include <ncode/net/net_common.h>
#include <stddef.h>
#include <random>
#include <set>
#include <string>
#include <tuple>

namespace ctr {
class PathProvider;
} /* namespace ctr */
namespace nc {
namespace lp {
class DemandMatrix;
} /* namespace lp */
namespace viz {
class HtmlPage;
} /* namespace viz */
} /* namespace nc */

namespace ctr {

// Each aggregate is identified by a combination of src, dst.
class AggregateId {
 public:
  AggregateId(nc::net::GraphNodeIndex src, nc::net::GraphNodeIndex dst)
      : src_(src), dst_(dst) {}

  explicit AggregateId(const std::pair<nc::net::GraphNodeIndex,
                                       nc::net::GraphNodeIndex>& src_and_dst)
      : src_(src_and_dst.first), dst_(src_and_dst.second) {}

  nc::net::GraphNodeIndex src() const { return src_; }

  nc::net::GraphNodeIndex dst() const { return dst_; }

  std::string ToString(const nc::net::GraphStorage& graph) const;

  // Returns the delay of the shortest path through the graph that this
  // aggregate can be routed on.
  nc::net::Delay GetSPDelay(const nc::net::GraphStorage& graph) const;

  // Returns the shortest path.
  std::unique_ptr<nc::net::Walk> GetSP(
      const nc::net::GraphStorage& graph) const;

  AggregateId Reverse() const { return {dst_, src_}; }

  friend bool operator<(const AggregateId& a, const AggregateId& b);
  friend bool operator==(const AggregateId& a, const AggregateId& b);
  friend bool operator!=(const AggregateId& a, const AggregateId& b);

 private:
  nc::net::GraphNodeIndex src_;
  nc::net::GraphNodeIndex dst_;
};

// A bandwidth demand and flow count.
using DemandAndFlowCount = std::pair<nc::net::Bandwidth, size_t>;

// A path and a fraction of demand that should go over it.
using RouteAndFraction = std::pair<const nc::net::Walk*, double>;

// For each aggregate, a bandwidth demand and a flow count.
class TrafficMatrix {
 public:
  // Generates a TM where each aggregate's flow count is proportional to its
  // traffic volume.
  static std::unique_ptr<TrafficMatrix> ProportionalFromDemandMatrix(
      const nc::lp::DemandMatrix& demand_matrix,
      size_t top_aggregate_flow_count = 1000);

  // Distributes a total number of flows based on traffic volume.
  static std::unique_ptr<TrafficMatrix> DistributeFromDemandMatrix(
      const nc::lp::DemandMatrix& demand_matrix,
      size_t total_flow_count = 10000);

  // Constructs an empty traffic matrix, or optionally let the caller
  // pre-populate it with demands and flow counts.
  explicit TrafficMatrix(
      const nc::net::GraphStorage* graph,
      const std::map<AggregateId, DemandAndFlowCount>& demand_and_counts = {})
      : graph_(graph) {
    demands_ = demand_and_counts;
  }

  // Constructs a traffic matrix from a demand matrix. The demand matrix has no
  // flow counts, so they need to be supplied explicitly per src/dst pair (the
  // second argument). If there is no entry in the map for a src/dst pair its
  // flow count is assumed to be 1.
  explicit TrafficMatrix(
      const nc::lp::DemandMatrix& demand_matrix,
      const std::map<nc::lp::SrcAndDst, size_t>& flow_counts = {});

  const std::map<AggregateId, DemandAndFlowCount>& demands() const {
    return demands_;
  }

  void AddDemand(const AggregateId& aggregate_id,
                 const DemandAndFlowCount& demand_and_flow_count);

  // Scales demands by a factor. If no aggregates are specified will scale all
  // aggregates.
  std::unique_ptr<TrafficMatrix> ScaleDemands(
      double factor, const std::set<AggregateId>& to_scale) const;

  // Adds the same value to all demands.
  std::unique_ptr<TrafficMatrix> AddToDemands(
      nc::net::Bandwidth value, const std::set<AggregateId>& to_extend) const;

  // A DemandMatrix is similar to TrafficMatrix, but has no flow counts.
  std::unique_ptr<nc::lp::DemandMatrix> ToDemandMatrix() const;

  // Returns a new traffic matrix with the same aggregates, but
  // 'aggregate_count' aggregates have their demand/flow count +- a fraction of
  // the demand/flow count in this one.
  std::unique_ptr<TrafficMatrix> Randomize(double demand_fraction,
                                           double flow_count_fraction,
                                           size_t aggregate_count,
                                           std::mt19937* rnd) const;

  const nc::net::GraphStorage* graph() const { return graph_; }

  // Dumps the entire TM to string.
  std::string ToString() const;

  // Dumps a single aggregate from the TM to string.
  std::string AggregateToString(const AggregateId& aggregate) const;

  // Prints a summary of the TM.
  std::string SummaryToString() const;

  std::pair<nc::net::Bandwidth, nc::net::Bandwidth> MinMaxAggregates() const;

 protected:
  // The graph.
  const nc::net::GraphStorage* graph_;

 private:
  // For each aggregate its demand and its flow count.
  std::map<AggregateId, DemandAndFlowCount> demands_;

  DISALLOW_COPY_AND_ASSIGN(TrafficMatrix);
};

// Fraction of flows that change from one path to another.
struct FlowPathChange {
  double fraction;
  const nc::net::Walk* from;
  const nc::net::Walk* to;
};

std::vector<FlowPathChange> GetFractionDeltas(
    const std::vector<RouteAndFraction>& prev,
    const std::vector<RouteAndFraction>& next);

// The difference between the same aggregate in two different outputs
// (RoutingConfigurations).
struct AggregateDelta {
  // List of changes.
  std::vector<FlowPathChange> changes;

  // New routes added (f == 0 -> f != 0).
  size_t routes_added = 0;

  // Old routes deleted (f != 0 -> f == 0).
  size_t routes_removed = 0;

  // Routes updated (f != 0 -> f != 0)
  size_t routes_updated = 0;

  double FractionDelta(double relative_change_threshold = 0.0) const;

  double FractionOnLongerPath(double relative_change_threshold = 0.0) const;
};

struct RoutingConfigurationDelta {
  // Fraction of total flows that changed path.
  double total_flow_fraction_delta;

  // Fraction of total volume that changed path.
  double total_volume_fraction_delta;

  // Fraction of total flows that changed to a longer path.
  double total_flow_fraction_on_longer_path;

  // Fraction of total volume that changed to a longer path.
  double total_volume_fraction_on_longer_path;

  // The sum of delays experienced by all flows is D. This is
  // (D_new - D_old) / D_old
  double total_per_flow_delay_delta;

  // Same as above, but absolute. D_new - D_old.
  nc::net::Delay total_per_flow_delay_delta_absolute;

  // Per-aggregate deltas.
  std::map<AggregateId, AggregateDelta> aggregates;

  // Returns the number of routed added, removed and updated
  std::tuple<size_t, size_t, size_t> TotalRoutes() const;

  std::string ToString(const nc::net::GraphStorage& graph) const;
};

// Extends a TM with for each aggregate a set of paths and a fraction of demand
// to route.
class RoutingConfiguration : public TrafficMatrix {
 public:
  static constexpr char kDefaultOptimizerString[] = "UNKNOWN";

  // Loads a routing configuration from the text produced by SerializeToText.
  static std::unique_ptr<RoutingConfiguration> LoadFromSerializedText(
      const TrafficMatrix& base_matrix,
      const std::vector<std::string>& node_order, const std::string& text,
      PathProvider* path_provider);

  explicit RoutingConfiguration(const TrafficMatrix& base_matrix)
      : TrafficMatrix(base_matrix.graph(), base_matrix.demands()),
        time_to_compute_(0),
        optimizer_string_(kDefaultOptimizerString) {}

  void AddRouteAndFraction(
      const AggregateId& aggregate_id,
      const std::vector<RouteAndFraction>& routes_and_fractions);

  const std::vector<RouteAndFraction>& FindRoutesOrDie(
      const AggregateId& aggregate_id) const;

  const std::map<AggregateId, std::vector<RouteAndFraction>>& routes() const {
    return configuration_;
  }

  std::string ToString() const;

  std::string AggregateToString(const AggregateId& aggregate,
                                double* aggregate_contribution = nullptr) const;

  // Renders the configuration as a graph on a web page.
  void ToHTML(nc::viz::HtmlPage* out) const;

  // A simple text-based serialization of the routing configuration. Each line
  // after the first one represents one path and has the following
  // comma-separated format:
  // node_0,...,node_n,fraction_of_demand
  // where node_0 is the index of the source node, node_n is the index of the
  // destination node (indices are based on the original .graph file) and
  // fraction_of_demand is the fraction of the total demand between the source
  // and the destination that goes along this path. This representation assumes
  // that there is at most one link between any two nodes in the graph. The
  // first line has the following format:
  // optimizer,compute_time
  // where num_paths is the number of paths, optimizer is a string that
  // identifies the optimizer that was used to generate this configuration and
  // compute_time is the time (in milliseconds) it took to produce it.
  std::string SerializeToText(const std::vector<std::string>& node_order) const;

  // Computes the difference between this routing configuration and another.
  // Both should have the same aggregates.
  RoutingConfigurationDelta GetDifference(
      const RoutingConfiguration& other) const;

  // For each aggregate, returns the indices of the aggregate's paths in the K
  // shortest paths sequence of the aggregate.
  std::map<AggregateId, std::vector<size_t>> GetKValues() const;

  // Returns sum of all flows' delay. If the sp argument is true will only add
  // up the flows' best possible shortest path delay.
  nc::net::Delay TotalPerFlowDelay(bool sp = false) const;

  // Returns a map from a link to its utilization.
  nc::net::GraphLinkMap<double> LinkUtilizations() const;

  // Number of aggregates that cross links with utilization > 1.
  size_t OverloadedAggregates() const;

  // Returns the max link utilization.
  double MaxLinkUtilization() const;

  // Returns the maximum number of paths an aggregate has.
  size_t MaxNumberOfPathsInAggregate() const {
    size_t out = 0;
    for (const auto& aggregate_and_routes : configuration_) {
      size_t count = aggregate_and_routes.second.size();
      out = std::max(out, count);
    }

    return out;
  }

  // Returns for each aggregate the number of paths in the aggregate.
  std::vector<size_t> NumberOfPathsInAggregate() const {
    std::vector<size_t> out;
    for (const auto& aggregate_and_routes : configuration_) {
      size_t count = aggregate_and_routes.second.size();
      out.emplace_back(count);
    }

    return out;
  }

  // Makes a copy.
  std::unique_ptr<RoutingConfiguration> Copy() const;

  void set_time_to_compute(std::chrono::milliseconds time) {
    time_to_compute_ = time;
  }

  void set_optimizer_string(const std::string& opt_string) {
    optimizer_string_ = opt_string;
  }

  const std::string& optimizer_string() const { return optimizer_string_; }

  const std::chrono::milliseconds& time_to_compute() const {
    return time_to_compute_;
  }

 private:
  // For each aggregate the path and fraction through the network.
  std::map<AggregateId, std::vector<RouteAndFraction>> configuration_;

  // Optional fields that identify the time it took to compute the solution and
  // the optimizer that was used. Set to 0 if unknown.
  std::chrono::milliseconds time_to_compute_;
  std::string optimizer_string_;
};

// A set of aggregates that share capacity.
struct AggregatesAndCapacity {
  // Aggregates and fractions.
  std::vector<std::pair<AggregateId, double>> aggregates;

  // Capacity.
  nc::net::Bandwidth capacity;
};

// For each aggregate, a set of other aggregates it shares capacity with.
class CompetingAggregates {
 public:
  void AddAggregatesAndCapacity(
      const AggregateId& aggregate_id,
      const std::vector<AggregatesAndCapacity>& aggregates_and_capacity);

  const std::map<AggregateId, std::vector<AggregatesAndCapacity>>& aggregates()
      const {
    return aggregates_;
  }

  void Clear() { aggregates_.clear(); }

 private:
  std::map<AggregateId, std::vector<AggregatesAndCapacity>> aggregates_;
};

// For an aggregate combines binned history and flow count. The history spans
// bin_size * bins.size() time.
class AggregateHistory {
 public:
  AggregateHistory(const std::vector<uint64_t>& bins,
                   std::chrono::milliseconds bin_size, uint64_t flow_count)
      : bins_(bins), bin_size_(bin_size), flow_count_(flow_count) {}

  // Generates a synthetic traffic matrix where all the bins have the same
  // level. The level is chosen so that it matches a given per-second rate.
  AggregateHistory(nc::net::Bandwidth rate, size_t bin_count,
                   std::chrono::milliseconds bin_size, size_t flow_count);

  AggregateHistory(const AggregateHistory& other)
      : AggregateHistory(other.bins_, other.bin_size_, other.flow_count_) {}

  // The mean rate.
  nc::net::Bandwidth mean_rate() const;

  // The max rate.
  nc::net::Bandwidth max_rate() const;

  // Maximum queue size at a given rate.
  std::chrono::milliseconds MaxQueueAtRate(nc::net::Bandwidth bandwidth) const;

  // Returns the max rate at which no bin will experience queue longer than
  // 'max_queue'.
  nc::net::Bandwidth MaxRateAtQueue(std::chrono::milliseconds max_queue) const;

  // Returns a vector with the per-second mean rates.
  std::vector<nc::net::Bandwidth> PerSecondMeans() const;

  // Adds the given rate to all bins.
  AggregateHistory AddRate(nc::net::Bandwidth rate) const;

  // Removes the given rate from all bins.
  AggregateHistory SubtractRate(nc::net::Bandwidth rate) const;

  const std::vector<uint64_t>& bins() const { return bins_; }

  void SetBin(size_t i, uint64_t value) { bins_[i] = value; }

  void set_flow_count(uint64_t flow_count) { flow_count_ = flow_count; }

  std::chrono::milliseconds bin_size() const { return bin_size_; }

  uint64_t flow_count() const { return flow_count_; }

 private:
  // Bytes transmitted per bin in this aggregate's history.
  std::vector<uint64_t> bins_;

  // How long each bin lasts.
  std::chrono::milliseconds bin_size_;

  // Number of flows.
  uint64_t flow_count_;
};

// Returns an aggregate history that will have a given mean rate.
ctr::AggregateHistory GetDummyHistory(nc::net::Bandwidth rate,
                                      std::chrono::milliseconds bin_size,
                                      std::chrono::milliseconds duration,
                                      size_t flow_count);

// Message from TLDR to the controller. This is essentially a request for
// capacity. For each aggregate that TLDR's endpoint controls there is the bins
// seen over the last period.
class TLDRRequest : public ::nc::htsim::Message {
 public:
  static constexpr uint8_t kTLDRRequestType = 181;

  TLDRRequest(nc::net::IPAddress ip_src, nc::net::IPAddress ip_dst,
              nc::EventQueueTime time_sent, uint64_t round_id, bool quick,
              const std::map<AggregateId, AggregateHistory>& aggregates);

  const std::map<AggregateId, AggregateHistory>& aggregates();

  nc::htsim::PacketPtr Duplicate() const override;

  std::string ToString() const override;

  uint64_t round_id() const { return round_id_; }

  bool quick() const { return quick_; }

 private:
  // Each request carries an id that identifies the round.
  uint64_t round_id_;

  // For each aggregate the history.
  std::map<AggregateId, AggregateHistory> aggregates_;

  // If this flag is set the optimization round will be quick -- only currently
  // installed paths in the network will be used. All requests for the same
  // round should carry the same quick_ flag.
  bool quick_;
};

// Message from the controller to TLDR, to force the generation of a
// TLDRRequest.
class TLDRForceRequest : public ::nc::htsim::Message {
 public:
  static constexpr uint8_t kTLDRForceRequestType = 170;

  TLDRForceRequest(nc::net::IPAddress ip_src, nc::net::IPAddress ip_dst,
                   nc::EventQueueTime time_sent);

  nc::htsim::PacketPtr Duplicate() const override;

  std::string ToString() const override;
};

// Sent from TLDR to the controller to trigger a re-optimization outside of the
// regular periodic optimizations. The controller will then generate
// TLDRForceRequests to all TLDR instances and they will send TLDRRequests.
class TLDRTriggerReoptimize : public ::nc::htsim::Message {
 public:
  static constexpr uint8_t kTLDRTriggerReoptimizeType = 171;

  TLDRTriggerReoptimize(nc::net::IPAddress ip_src, nc::net::IPAddress ip_dst,
                        nc::EventQueueTime time_sent);

  nc::htsim::PacketPtr Duplicate() const override;

  std::string ToString() const override;
};

struct PathUpdateState {
  PathUpdateState(nc::htsim::PacketTag tag,
                  nc::net::DevicePortNumber output_port_num,
                  const RouteAndFraction& route_and_fraction)
      : tag(tag),
        output_port_num(output_port_num),
        route_and_fraction(route_and_fraction) {}

  // Packets that belong to this path are identified by the tag.
  nc::htsim::PacketTag tag;

  // Port on the switch managed by the observer that this path's traffic should
  // exit via.
  nc::net::DevicePortNumber output_port_num;

  // The route and the fraction of the aggregate that should go on it.
  RouteAndFraction route_and_fraction;
};

// Aggregate state sent from the controller to the observer.
struct AggregateUpdateState {
  AggregateUpdateState(AggregateId aggregate_id, ::nc::htsim::MatchRuleKey key)
      : aggregate_id(aggregate_id), key(key) {}

  AggregateId aggregate_id;

  // Incoming traffic should match this key to be considered part of the
  // aggregate.
  ::nc::htsim::MatchRuleKey key;

  // Paths in the aggregate.
  std::vector<PathUpdateState> paths;

  // For each path in paths, a set of aggregates the path competes with and
  // capacity it competed for. Used when determining whether to trigger
  // optimization or not. If empty will never trigger optimizations.
  std::vector<AggregatesAndCapacity> competing_aggregates;
};

// Message sent from the controller to the observer. Has for each aggregate the
// paths and rates to enforce,
class TLDRUpdate : public ::nc::htsim::Message {
 public:
  static constexpr uint8_t kTLDRUpdateType = 182;

  TLDRUpdate(
      nc::net::IPAddress ip_src, nc::net::IPAddress ip_dst,
      nc::EventQueueTime time_sent,
      const std::map<AggregateId, AggregateUpdateState>& aggregate_updates,
      const std::map<AggregateId, AggregateHistory>& additional_histories);

  const std::map<AggregateId, AggregateUpdateState>& aggregates() const {
    return aggregates_;
  }

  const std::map<AggregateId, AggregateHistory>& additional_histories() const {
    return additional_histories_;
  }

  nc::htsim::PacketPtr Duplicate() const override;

  std::string ToString() const override;

 private:
  // Each aggregate's paths. There will be an AggregateUpdateState for each
  // aggregate controlled by the TLDR this update is for.
  std::map<AggregateId, AggregateUpdateState> aggregates_;

  // Histories of aggregates that compete with any of the aggregates_'s paths.
  // Used when determining when to trigger updates.
  std::map<AggregateId, AggregateHistory> additional_histories_;
};

}  // namespace ctr
#endif
