#include "metrics_test_util.h"

#include <cstdlib>
#include <iostream>
#include <string>

#include "ncode/common.h"

namespace nc {
namespace metrics {
namespace test {

MetricFixtureBase::MetricFixtureBase() {
  metric_manager_ = make_unique<MetricManager>();
  metric_manager_->SetOutput(kTestOutput);
  timestamp_provider_ = metric_manager_->timestamp_provider();
}

void MetricFixtureBase::TearDownBase() {
  std::string cmd(std::string("rm -rf ") + kTestOutput);
  auto res = system(cmd.c_str());
  Unused(res);
}

}  // namespace metrics
}  // namespace test
}  // namespace ncode
