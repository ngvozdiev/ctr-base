#include "opt.h"

#include <gtest/gtest.h>
#include <set>

namespace ctr {

constexpr nc::net::Bandwidth kDefaultLinkspeed =
    nc::net::Bandwidth::FromBitsPerSecond(1000000);

static nc::net::GraphBuilder GetGraph() {
  using namespace std::chrono;
  using namespace nc::net;
  GraphBuilder builder;
  builder.AddLink({"C", "D", kDefaultLinkspeed, milliseconds(10)});
  builder.AddLink({"D", "C", kDefaultLinkspeed, milliseconds(10)});
  builder.AddLink({"A", "B", kDefaultLinkspeed, milliseconds(10)});
  builder.AddLink({"B", "A", kDefaultLinkspeed, milliseconds(10)});
  builder.AddLink({"A", "D", kDefaultLinkspeed, milliseconds(6)});
  builder.AddLink({"D", "A", kDefaultLinkspeed, milliseconds(6)});
  builder.AddLink({"B", "C", kDefaultLinkspeed, milliseconds(5)});
  builder.AddLink({"C", "B", kDefaultLinkspeed, milliseconds(5)});

  return builder;
}

class TestBase : public ::testing::Test {
 protected:
  TestBase() : graph_(GetGraph()) {
    path_provider_ = nc::make_unique<PathProvider>(&graph_);
  }

  // Given an output and a path will check if the output contains the path and
  // sends a given fraction down that path.
  bool HasPath(const RoutingConfiguration& output, const std::string& path_str,
               double fraction) {
    std::unique_ptr<nc::net::Walk> path = graph_.WalkFromStringOrDie(path_str);
    AggregateId id(path->FirstHop(graph_), path->LastHop(graph_));

    for (const auto& aggregate_and_routes : output.routes()) {
      const AggregateId& aggregate_id = aggregate_and_routes.first;
      if (aggregate_id != id) {
        continue;
      }

      const std::vector<RouteAndFraction>& routes = aggregate_and_routes.second;
      for (const auto& route_and_fraction : routes) {
        if (*route_and_fraction.first == *path) {
          double fraction_in_result = route_and_fraction.second;
          if (std::fabs(fraction_in_result - fraction) < 0.01) {
            return true;
          }
        }
      }
    }

    return false;
  }

  void AddAggregate(const std::string& src, const std::string& dst,
                    uint32_t num_flows, nc::net::Bandwidth total_volume) {
    AggregateId id(graph_.NodeFromStringOrDie(src),
                   graph_.NodeFromStringOrDie(dst));
    tm_.AddDemand(id, {total_volume, num_flows});
  }

  std::unique_ptr<PathProvider> path_provider_;
  nc::net::GraphStorage graph_;
  TrafficMatrix tm_;
};

class ShortestPathTest : public TestBase {
 protected:
  ShortestPathTest() : sp_optimizer_(std::move(path_provider_)) {}

  void AddSPAggregate(const std::string& src, const std::string& dst) {
    AddAggregate(src, dst, 1, kDefaultLinkspeed);
  }

  ShortestPathOptimizer sp_optimizer_;
};

TEST_F(ShortestPathTest, SingleAggregate) {
  AddSPAggregate("A", "B");
  auto routing = sp_optimizer_.Optimize(tm_);
  ASSERT_TRUE(HasPath(*routing, "[A->B]", 1.0));
}

TEST_F(ShortestPathTest, TwoAggregates) {
  ShortestPathOptimizer sp_optimizer(std::move(path_provider_));

  AddSPAggregate("A", "B");
  AddSPAggregate("A", "C");
  auto routing = sp_optimizer_.Optimize(tm_);
  ASSERT_TRUE(HasPath(*routing, "[A->B]", 1.0));
  ASSERT_TRUE(HasPath(*routing, "[A->B, B->C]", 1.0));
}

class B4Test : public TestBase {
 protected:
  B4Test() : b4_optimizer_(std::move(path_provider_)) {}

  void AddB4Aggregate(const std::string& src, const std::string& dst,
                      nc::net::Bandwidth bw) {
    AddAggregate(src, dst, 1, bw);
  }

  B4Optimizer b4_optimizer_;
};

// Fits on the shortest path.
TEST_F(B4Test, SinglePath) {
  AddB4Aggregate("A", "B", kDefaultLinkspeed.bps());

  Input input(1.0, std::move(input_map_));
  B4Optimizer optimizer(&path_cache_);
  B4Output output = optimizer.Optimize(input);
  ASSERT_FALSE(output.IsOversubscribed());

  ASSERT_EQ(1ul, output.aggregates().size());
  const AggregateOutput& aggregate_output =
      ncode::FindOrDie(output.aggregates(), 10);

  // +1 for the backup path.
  ASSERT_EQ(2ul, aggregate_output.paths().size());

  const PathOutput* path_output = aggregate_output.PathsSorted().front();
  ASSERT_EQ(1.0, path_output->fraction());

  const ncode::net::GraphPath* graph_path =
      graph_storage_.PathFromStringOrDie("[A->B]", 10);
  ASSERT_EQ(graph_path, path_output->path());
}

// Does not fit on the shortest path
TEST_F(B4Test, SinglePathNoFit) {
  AddAggregate("A", "B", 10, 100, kDefaultLinkspeed.bps() * 1.2);

  Input input(0.9, std::move(input_map_));
  B4Optimizer optimizer(&path_cache_);
  B4Output output = optimizer.Optimize(input);
  ASSERT_FALSE(output.IsOversubscribed());

  ASSERT_EQ(1ul, output.aggregates().size());
  const AggregateOutput& aggregate_output =
      ncode::FindOrDie(output.aggregates(), 10);
  ASSERT_EQ(2ul, aggregate_output.paths().size());

  // This should be the shorter of the 2 paths (A->B).
  const PathOutput* path_output = aggregate_output.PathsSorted().front();
  ASSERT_DOUBLE_EQ(0.9, path_output->fraction() * 1.2);

  const ncode::net::GraphPath* graph_path =
      graph_storage_.PathFromStringOrDie("[A->B]", 10);
  ASSERT_EQ(graph_path, path_output->path());

  const PathOutput* other_path_output = aggregate_output.PathsSorted().back();
  ASSERT_DOUBLE_EQ(0.3, other_path_output->fraction() * 1.2);

  graph_path = graph_storage_.PathFromStringOrDie("[A->D, D->C, C->B]", 10);
  ASSERT_EQ(graph_path, other_path_output->path());
}

// A couple of aggregates that do not fit on the shortest path
TEST_F(B4Test, MultiSinglePathNoFit) {
  AddAggregate("A", "B", 10, 100, kDefaultLinkspeed.bps() * 0.6);
  AddAggregate("A", "C", 11, 100, kDefaultLinkspeed.bps() * 0.6);

  Input input(0.9, std::move(input_map_));
  B4Optimizer optimizer(&path_cache_);
  B4Output output = optimizer.Optimize(input);
  ASSERT_FALSE(output.IsOversubscribed());
  ASSERT_EQ(2ul, output.aggregates().size());

  ASSERT_TRUE(HasPath(output, "[A->B]", 0.9 / 1.2, 10));
  ASSERT_TRUE(HasPath(output, "[A->D, D->C, C->B]", 0.3 / 1.2, 10));
  ASSERT_TRUE(HasPath(output, "[A->B, B->C]", 0.9 / 1.2, 11));
  ASSERT_TRUE(HasPath(output, "[A->D, D->C]", 0.3 / 1.2, 11));
}

// A couple of aggregates, one does not fit.
TEST_F(B4Test, MultiSinglePathOneNoFit) {
  AddAggregate("A", "B", 10, 100, kDefaultLinkspeed.bps() * 0.3);
  AddAggregate("A", "C", 11, 100, kDefaultLinkspeed.bps() * 0.9);

  Input input(0.9, std::move(input_map_));
  B4Optimizer optimizer(&path_cache_);
  B4Output output = optimizer.Optimize(input);
  ASSERT_FALSE(output.IsOversubscribed());
  ASSERT_EQ(2ul, output.aggregates().size());

  // Should be the same as in the previous test!
  ASSERT_TRUE(HasPath(output, "[A->B]", 0.9 / 1.2, 10));
  ASSERT_TRUE(HasPath(output, "[A->D, D->C, C->B]", 0.3 / 1.2, 10));
  ASSERT_TRUE(HasPath(output, "[A->B, B->C]", 0.9 / 1.2, 11));
  ASSERT_TRUE(HasPath(output, "[A->D, D->C]", 0.3 / 1.2, 11));
}

// Traffic does not fit.
TEST_F(B4Test, NoFit) {
  AddAggregate("A", "B", 10, 100, kDefaultLinkspeed.bps() * 3);

  Input input(1.0, std::move(input_map_));
  B4Optimizer optimizer(&path_cache_);
  B4Output output = optimizer.Optimize(input);
  ASSERT_TRUE(output.IsOversubscribed());
}

}  // namespace ctr
