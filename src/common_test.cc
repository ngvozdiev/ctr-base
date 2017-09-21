#include "common.h"

#include <gtest/gtest.h>

namespace ctr {

class FractionDeltaFixture : public ::testing::Test {
 protected:
  const nc::net::Walk* GetPath() {
    std::unique_ptr<nc::net::Walk> new_path = nc::make_unique<nc::net::Walk>();
    paths_.emplace_back(std::move(new_path));
    return paths_.back().get();
  }

  std::vector<std::unique_ptr<nc::net::Walk>> paths_;
};

TEST_F(FractionDeltaFixture, Simple) {
  std::vector<RouteAndFraction> from;
  std::vector<RouteAndFraction> to;
  ASSERT_TRUE(GetFractionDeltas(from, to).empty());
}

TEST_F(FractionDeltaFixture, NoChangeSingle) {
  auto p1 = GetPath();

  std::vector<RouteAndFraction> from = {{p1, 1.0}};
  std::vector<RouteAndFraction> to = {{p1, 1.0}};

  std::vector<FlowPathChange> changes = GetFractionDeltas(from, to);
  ASSERT_EQ(1ul, changes.size());

  const FlowPathChange& change_one = changes.front();
  ASSERT_EQ(p1, change_one.from);
  ASSERT_EQ(p1, change_one.to);
  ASSERT_DOUBLE_EQ(1.0, change_one.fraction);
}

TEST_F(FractionDeltaFixture, NoChangeMulti) {
  auto p1 = GetPath();
  auto p2 = GetPath();

  std::vector<RouteAndFraction> from = {{p1, 0.2}, {p2, 0.8}};
  std::vector<RouteAndFraction> to = {{p1, 0.2}, {p2, 0.8}};

  std::vector<FlowPathChange> changes = GetFractionDeltas(from, to);
  ASSERT_EQ(2ul, changes.size());

  const FlowPathChange& change_one = changes.front();
  ASSERT_EQ(p1, change_one.from);
  ASSERT_EQ(p1, change_one.to);
  ASSERT_DOUBLE_EQ(0.2, change_one.fraction);

  const FlowPathChange& change_two = changes.back();
  ASSERT_EQ(p2, change_two.from);
  ASSERT_EQ(p2, change_two.to);
  ASSERT_DOUBLE_EQ(0.8, change_two.fraction);
}

TEST_F(FractionDeltaFixture, Flip) {
  auto p1 = GetPath();
  auto p2 = GetPath();

  std::vector<RouteAndFraction> from = {{p1, 0.2}, {p2, 0.8}};
  std::vector<RouteAndFraction> to = {{p1, 0.8}, {p2, 0.2}};

  std::vector<FlowPathChange> changes = GetFractionDeltas(from, to);
  ASSERT_EQ(3ul, changes.size());

  const FlowPathChange& change_one = changes.front();
  ASSERT_EQ(p1, change_one.from);
  ASSERT_EQ(p1, change_one.to);
  ASSERT_DOUBLE_EQ(0.2, change_one.fraction);

  const FlowPathChange& change_two = changes[1];
  ASSERT_EQ(p2, change_two.from);
  ASSERT_EQ(p1, change_two.to);
  ASSERT_DOUBLE_EQ(0.6, change_two.fraction);

  const FlowPathChange& change_three = changes.back();
  ASSERT_EQ(p2, change_three.from);
  ASSERT_EQ(p2, change_three.to);
  ASSERT_DOUBLE_EQ(0.2, change_three.fraction);
}

TEST_F(FractionDeltaFixture, Changes) {
  auto p1 = GetPath();
  auto p2 = GetPath();
  auto p3 = GetPath();

  std::vector<RouteAndFraction> from = {{p1, 1.0}};
  std::vector<RouteAndFraction> to = {{p2, 0.7}, {p3, 0.3}};

  std::vector<FlowPathChange> changes = GetFractionDeltas(from, to);
  ASSERT_EQ(2ul, changes.size());

  const FlowPathChange& change_one = changes.front();
  ASSERT_EQ(p1, change_one.from);
  ASSERT_EQ(p2, change_one.to);
  ASSERT_DOUBLE_EQ(0.7, change_one.fraction);

  const FlowPathChange& change_two = changes.back();
  ASSERT_EQ(p1, change_two.from);
  ASSERT_EQ(p3, change_two.to);
  ASSERT_DOUBLE_EQ(0.3, change_two.fraction);
}

TEST_F(FractionDeltaFixture, ChangesTwo) {
  auto p1 = GetPath();
  auto p2 = GetPath();
  auto p3 = GetPath();

  std::vector<RouteAndFraction> from = {{p2, 0.7}, {p3, 0.3}};
  std::vector<RouteAndFraction> to = {{p1, 1.0}};

  std::vector<FlowPathChange> changes = GetFractionDeltas(from, to);
  ASSERT_EQ(2ul, changes.size());

  const FlowPathChange& change_one = changes.front();
  ASSERT_EQ(p2, change_one.from);
  ASSERT_EQ(p1, change_one.to);
  ASSERT_DOUBLE_EQ(0.7, change_one.fraction);

  const FlowPathChange& change_two = changes.back();
  ASSERT_EQ(p3, change_two.from);
  ASSERT_EQ(p1, change_two.to);
  ASSERT_DOUBLE_EQ(0.3, change_two.fraction);
}

TEST_F(FractionDeltaFixture, NewPath) {
  auto p1 = GetPath();
  auto p2 = GetPath();
  auto p3 = GetPath();
  auto p4 = GetPath();

  std::vector<RouteAndFraction> from = {{p4, 0.1}, {p1, 0.7}, {p3, 0.2}};
  std::vector<RouteAndFraction> to = {{p4, 0.1}, {p1, 0.6}, {p2, 0.3}};

  std::vector<FlowPathChange> changes = GetFractionDeltas(from, to);
  ASSERT_EQ(4ul, changes.size());

  const FlowPathChange& change_one = changes.front();
  ASSERT_EQ(p4, change_one.from);
  ASSERT_EQ(p4, change_one.to);
  ASSERT_DOUBLE_EQ(0.1, change_one.fraction);

  const FlowPathChange& change_two = changes[1];
  ASSERT_EQ(p1, change_two.from);
  ASSERT_EQ(p1, change_two.to);
  ASSERT_DOUBLE_EQ(0.6, change_two.fraction);

  const FlowPathChange& change_three = changes[2];
  ASSERT_EQ(p1, change_three.from);
  ASSERT_EQ(p2, change_three.to);
  ASSERT_DOUBLE_EQ(0.1, change_three.fraction);

  const FlowPathChange& change_four = changes.back();
  ASSERT_EQ(p3, change_four.from);
  ASSERT_EQ(p2, change_four.to);
  ASSERT_DOUBLE_EQ(0.2, change_four.fraction);
}

}  // namespace ctr
