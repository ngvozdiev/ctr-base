#include "dist_model.h"

#include "gmock/gmock.h"
#include <algorithm>
#include <chrono>
#include <memory>

#include "../common.h"

namespace ctr {
namespace {

using ::testing::ElementsAre;
using ::testing::Pair;
using ::testing::DoubleEq;
using ::testing::Eq;

static size_t kDefaultInitialQuantization = 10;
static std::chrono::milliseconds kDefaultBinSize =
    std::chrono::milliseconds(100);
static size_t kDefaultFlowCount = 100;

class ProbTest : public ::testing::Test {
 protected:
  ProbTest() : distribution_({1001, 1002, 1003, 1004, 1005, 1006}) {}

  Distribution distribution_;
};

TEST_F(ProbTest, SimpleProbs) {
  ASSERT_THAT(
      distribution_.GetProbabilities(),
      ElementsAre(Pair(1001, DoubleEq(1 / 6.0)), Pair(1002, DoubleEq(1 / 6.0)),
                  Pair(1003, DoubleEq(1 / 6.0)), Pair(1004, DoubleEq(1 / 6.0)),
                  Pair(1005, DoubleEq(1 / 6.0)),
                  Pair(1006, DoubleEq(1 / 6.0))));
}

TEST_F(ProbTest, Percentiles) {
  ASSERT_THAT(distribution_.Percentile(1.0), DoubleEq(1006));
  ASSERT_THAT(distribution_.Percentile(0), DoubleEq(1001));
  ASSERT_THAT(distribution_.Percentile(0.5), DoubleEq(1003));
}

TEST_F(ProbTest, DataSeries) {
  ASSERT_THAT(
      distribution_.GetDataSeries(true),
      ElementsAre(Pair(1001, DoubleEq(1 / 6.0)), Pair(1002, DoubleEq(2 / 6.0)),
                  Pair(1003, DoubleEq(3 / 6.0)), Pair(1004, DoubleEq(4 / 6.0)),
                  Pair(1005, DoubleEq(5 / 6.0)), Pair(1006, DoubleEq(1))));
}

static ProbModelConfig GetDefaultConfig() {
  ProbModelConfig config;
  config.initial_quantization = kDefaultInitialQuantization;
  config.distribution_levels = 1024;
  return config;
}

class ProbModelFixture : public ::testing::Test {
 protected:
  ProbModelFixture() : model_(GetDefaultConfig()) {}

  void AddAggregate(uint64_t id, const std::vector<uint64_t>& bins_bytes) {
    auto new_history = nc::make_unique<AggregateHistory>(
        bins_bytes, kDefaultBinSize, kDefaultFlowCount);
    model_.AddAggregate(id, new_history.get());
    histories_.emplace_back(std::move(new_history));
  }

  // A number of uniformly distributed bins in a range.
  std::vector<uint64_t> GenerateBins(uint64_t from, uint64_t to,
                                     uint64_t count) {
    std::vector<uint64_t> out;
    auto dist = std::uniform_int_distribution<uint64_t>(from, to);
    for (size_t i = 0; i < count; ++i) {
      out.emplace_back(dist(rnd_));
    }

    return out;
  }

  bool OutputOk(const std::vector<bool>& model,
                const std::vector<ProbModelReply>& replies) {
    if (replies.size() != model.size()) {
      return false;
    }
    for (size_t i = 0; i < replies.size(); ++i) {
      if (replies[i].fit != model[i]) {
        return false;
      }
    }

    return true;
  }

  std::mt19937 rnd_;
  ProbModel model_;
  std::vector<std::unique_ptr<AggregateHistory>> histories_;
};

TEST_F(ProbModelFixture, Empty) { ASSERT_TRUE(model_.Query({}).empty()); }

TEST_F(ProbModelFixture, SingleAggregate) {
  AddAggregate(1, GenerateBins(1000, 100000, 600));

  nc::net::Bandwidth rate =
      nc::net::Bandwidth::FromBitsPerSecond(100000 * 10 * 8);
  nc::net::Bandwidth rate_too_low =
      nc::net::Bandwidth::FromBitsPerSecond(100000 * 8);

  std::vector<ProbModelReply> out =
      model_.Query({{ProbModelQuery::BOTH, {{1, 1.0}}, rate},
                    {ProbModelQuery::BOTH, {{1, 1.0}}, rate_too_low}});
  std::vector<bool> model = {true, false};
  ASSERT_TRUE(OutputOk(model, out));
}

TEST_F(ProbModelFixture, SingleAggregateSplit) {
  AddAggregate(1, GenerateBins(1000, 100000, 600));
  nc::net::Bandwidth rate =
      nc::net::Bandwidth::FromBitsPerSecond(100000 * 10 * 4);

  std::vector<ProbModelReply> out =
      model_.Query({{ProbModelQuery::BOTH, {{1, 0.5}}, rate}});
  std::vector<bool> model = {true};
  ASSERT_TRUE(OutputOk(model, out));
}

TEST_F(ProbModelFixture, SingleAggregateSplitTimings) {
  using namespace std::chrono;

  AddAggregate(1, GenerateBins(1000, 100000, 600));
  nc::net::Bandwidth rate =
      nc::net::Bandwidth::FromBitsPerSecond(100000 * 10 * 4);

  std::vector<ProbModelTimingData> timing_data;
  model_.Query({{ProbModelQuery::BOTH, {{1, 0.5}}, rate}}, &timing_data);
  ASSERT_EQ(1ul, timing_data.size());
  ASSERT_NE(microseconds::zero(), timing_data[0].queue_simulation);
  ASSERT_NE(microseconds::zero(), timing_data[0].split_aggregates);
  ASSERT_NE(microseconds::zero(), timing_data[0].total);
  ASSERT_NE(microseconds::zero(), timing_data[0].convolution_timing.binning);
  ASSERT_NE(microseconds::zero(), timing_data[0].convolution_timing.fft);
  ASSERT_NE(microseconds::zero(), timing_data[0].convolution_timing.ifft);
}

TEST_F(ProbModelFixture, MultiAggregate) {
  // One hundred random aggregates.
  for (size_t i = 0; i < 10000; ++i) {
    AddAggregate(i, GenerateBins(1000, 100000, 600));
  }

  uint64_t mean_bytes_per_bin = (1000 + 100000) / 2.0;

  nc::net::Bandwidth single_rate =
      nc::net::Bandwidth::FromBitsPerSecond(100000ul * 10 * 8);
  nc::net::Bandwidth large_rate = nc::net::Bandwidth::FromBitsPerSecond(
      (mean_bytes_per_bin * 1.1) * 10 * 8 * 10000);
  nc::net::Bandwidth small_rate =
      nc::net::Bandwidth::FromBitsPerSecond(1000 * 10 * 8 * 10000ul);

  std::vector<ProbModelReply> out =
      model_.Query({{ProbModelQuery::BOTH, {{1, 1.0}}, single_rate}});
  std::vector<bool> model = {true};
  ASSERT_TRUE(OutputOk(model, out));

  std::vector<std::pair<uint64_t, double>> all_aggregates;
  for (size_t i = 0; i < 10000; ++i) {
    all_aggregates.emplace_back(i, 1.0);
  }

  std::vector<ProbModelTimingData> timing_data;
  out = model_.Query({{ProbModelQuery::BOTH, all_aggregates, small_rate},
                      {ProbModelQuery::BOTH, all_aggregates, large_rate}},
                     &timing_data);

  for (size_t i = 0; i < out.size(); ++i) {
    LOG(INFO) << timing_data[i].ToString();
  }

  model = {false, true};
  ASSERT_TRUE(OutputOk(model, out));
}

}  // namespace
}  // namespace ctr
