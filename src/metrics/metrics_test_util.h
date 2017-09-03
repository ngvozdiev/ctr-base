#ifndef NCODE_METRICS_METRIC_TEST_UTIL_H
#define NCODE_METRICS_METRIC_TEST_UTIL_H

#include <cstdint>
#include <memory>

#include "gtest/gtest.h"
#include "metrics.h"

namespace nc {
namespace metrics {
namespace test {

static constexpr char kMetricId[] = "some_metric";
static constexpr char kMetricDesc[] = "Some description";
static constexpr char kMetricFieldOneDesc[] = "a field";
static constexpr char kMetricFieldTwoDesc[] = "another field";
static constexpr char kMetricFieldStrValue[] = "field str value";
static constexpr char kMetricAnotherFieldStrValue[] = "another field str value";
static constexpr char kTestOutput[] = "test_output";
static constexpr char kJunk[] = "junk";
static constexpr uint64_t kMetricFieldIntValue = 45;

class MetricFixtureBase {
 protected:
  MetricFixtureBase();

  void TearDownBase();

  const TimestampProviderInterface* timestamp_provider_;
  std::unique_ptr<MetricManager> metric_manager_;
};

// A MetricFixture that sets up a single output file.
class MetricFixture : public MetricFixtureBase, public ::testing::Test {
 protected:
  MetricFixture() : MetricFixtureBase() {}
  void TearDown() override { TearDownBase(); }
};

}  // namespace test
}  // namespace metrics
}  // namespace ncode

#endif /* METRICS_METRIC_TEST_UTIL_H */
