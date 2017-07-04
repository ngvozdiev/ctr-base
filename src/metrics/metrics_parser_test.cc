#include "metrics_parser.h"

#include <stddef.h>
#include <numeric>

#include "gtest/gtest.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/substitute.h"
#include "metrics_test_util.h"

using namespace nc::metrics::parser;

namespace nc {
namespace metrics {
namespace test {

using DoubleProcessor = QueryCallbackProcessor<double, PBManifestEntry::DOUBLE>;

TEST_F(MetricFixture, DoubleParser) {
  auto* metric = metric_manager_->GetUnsafeMetric<double, std::string>(
      kMetricComonentId, kMetricDesc, kMetricFieldOneDesc);
  auto* handle = metric->GetHandle(kMetricFieldStrValue);
  for (size_t i = 0; i < 1000000; ++i) {
    handle->AddValue(i);
  }
  metric_manager_.reset();

  std::vector<double> all_values;
  DoubleProcessor::Callback double_callback = [&all_values](
      const Entry<double>& entry, const PBManifestEntry& manifest_entry,
      uint32_t manifest_index) {
    ASSERT_EQ(0ul, manifest_index);
    ASSERT_EQ(kMetricComonentId, manifest_entry.id());
    all_values.emplace_back(entry.value);
  };

  auto double_processor =
      make_unique<DoubleProcessor>(".*", kMetricFieldStrValue, double_callback);

  MetricsParser parser(kTestOutput);
  parser.AddProcessor(std::move(double_processor));
  parser.Parse();

  ASSERT_EQ(1000000ul, all_values.size());
  for (size_t i = 0; i < 1000000; ++i) {
    ASSERT_EQ(all_values[i], i);
  }
}

TEST_F(MetricFixture, ExternalNoMetrics) {
  auto* metric = metric_manager_->GetUnsafeMetric<uint64_t, std::string>(
      kMetricComonentId, kMetricDesc, kMetricFieldOneDesc);
  auto* handle = metric->GetHandle(kMetricFieldStrValue);
  for (size_t i = 0; i < 1000000; ++i) {
    handle->AddValue(i);
  }
  metric_manager_.reset();

  NumericMetricsResultHandle* result_handle = MetricsParserParse(
      kTestOutput, ".*", kJunk, 0, std::numeric_limits<uint64_t>::max(), 0);
  ASSERT_NE(nullptr, result_handle);
  ASSERT_FALSE(result_handle->Advance());
  MetricsParserResultHandleFree(result_handle);

  result_handle =
      MetricsParserParse(kTestOutput, kJunk, kMetricFieldStrValue, 0,
                         std::numeric_limits<uint64_t>::max(), 0);
  ASSERT_NE(nullptr, result_handle);
  ASSERT_FALSE(result_handle->Advance());
  MetricsParserResultHandleFree(result_handle);
}

TEST_F(MetricFixture, External) {
  auto* metric = metric_manager_->GetUnsafeMetric<uint64_t, std::string>(
      kMetricComonentId, kMetricDesc, kMetricFieldOneDesc);
  auto* handle = metric->GetHandle(kMetricFieldStrValue);
  for (size_t i = 0; i < 1000000; ++i) {
    handle->AddValue(i);
  }
  metric_manager_.reset();

  NumericMetricsResultHandle* result_handle =
      MetricsParserParse(kTestOutput, kMetricComonentId, kMetricFieldStrValue,
                         0, std::numeric_limits<uint64_t>::max(), 0);
  ASSERT_NE(nullptr, result_handle);
  ASSERT_TRUE(result_handle->Advance());

  size_t size = MetricsParserResultHandleSize(result_handle);
  ASSERT_EQ(1000000ul, size);

  std::vector<double> values_out(size);
  std::vector<uint64_t> timestamps_out(size);

  MetricsParserResultHandleCopyInto(result_handle, &timestamps_out[0],
                                    &values_out[0]);
  for (size_t i = 0; i < 1000000; ++i) {
    ASSERT_EQ(values_out[i], i);
    ASSERT_GT(timestamps_out[i], 0ul);
  }

  ASSERT_FALSE(result_handle->Advance());
  MetricsParserResultHandleFree(result_handle);
}

TEST_F(MetricFixture, BytesManifestSummary) {
  auto* metric =
      metric_manager_->GetUnsafeMetric<BytesBlob, std::string, uint64_t>(
          kMetricComonentId, kMetricDesc, kMetricFieldOneDesc,
          kMetricFieldTwoDesc);
  auto* handle = metric->GetHandle(kMetricFieldStrValue, 0);

  std::string v1 = "AA";
  for (size_t i = 0; i < 1000000; ++i) {
    BytesBlob bytes_blob;
    bytes_blob.set_bytes_value(v1);
    handle->AddValue(bytes_blob);
  }
  metric_manager_.reset();

  MetricsParser parser(kTestOutput);
  std::string manifest_summary = parser.ParseManifest().FullToString();

  ASSERT_NE(std::string::npos, manifest_summary.find(kMetricComonentId));
  ASSERT_NE(std::string::npos, manifest_summary.find(kMetricFieldOneDesc));
  ASSERT_NE(std::string::npos, manifest_summary.find(std::to_string(1000000)));
}

TEST_F(MetricFixture, ExternalBytes) {
  auto* metric = metric_manager_->GetUnsafeMetric<BytesBlob, std::string>(
      kMetricComonentId, kMetricDesc, kMetricFieldOneDesc);
  auto* handle = metric->GetHandle(kMetricFieldStrValue);

  std::string v1 = "AA";
  std::string v2 = "BBBB";

  for (size_t i = 0; i < 1000000; ++i) {
    BytesBlob bytes_blob;
    if (i % 2 == 0) {
      bytes_blob.set_bytes_value(v1);
    } else {
      bytes_blob.set_bytes_value(v2);
    }
    handle->AddValue(bytes_blob);
  }
  metric_manager_.reset();

  NumericMetricsResultHandle* numeric_result_handle =
      MetricsParserParse(kTestOutput, ".*", kMetricFieldStrValue, 0,
                         std::numeric_limits<uint64_t>::max(), 0);
  ASSERT_FALSE(numeric_result_handle->Advance());

  BytesMetricsResultHandle* bytes_result_handle =
      MetricsParserBytesParse(kTestOutput, ".*", kMetricFieldStrValue, 0,
                              std::numeric_limits<uint64_t>::max(), 0);
  ASSERT_NE(nullptr, bytes_result_handle);
  ASSERT_TRUE(bytes_result_handle->Advance());

  size_t size = MetricsParserBytesResultHandleSize(bytes_result_handle);
  ASSERT_EQ(1000000ul, size);

  size_t total_size = 0;
  for (size_t i = 0; i < 1000000; ++i) {
    uint64_t buffer_size = bytes_result_handle->BufferSize(i);
    if (i % 2 == 0) {
      ASSERT_EQ(2ul, buffer_size);
    } else {
      ASSERT_EQ(4ul, buffer_size);
    }

    total_size += buffer_size;
  }

  auto buffer = new char[total_size];
  std::vector<uint64_t> timestamps_out(size);
  MetricsParserBytesResultHandleCopyInto(bytes_result_handle,
                                         &timestamps_out[0], buffer);

  // At this point the buffer should contain all 1M objects back to back.
  size_t offset = 0;
  for (size_t i = 0; i < 1000000; ++i) {
    if (i % 2 == 0) {
      ASSERT_EQ('A', buffer[offset + 0]);
      ASSERT_EQ('A', buffer[offset + 1]);
      offset += 2;
    } else {
      ASSERT_EQ('B', buffer[offset + 0]);
      ASSERT_EQ('B', buffer[offset + 1]);
      ASSERT_EQ('B', buffer[offset + 2]);
      ASSERT_EQ('B', buffer[offset + 3]);
      offset += 4;
    }
  }

  MetricsParserBytesResultHandleFree(bytes_result_handle);
  delete[] buffer;
}

TEST_F(MetricFixture, ExternalTimestampLimits) {
  auto* metric = metric_manager_->GetUnsafeMetric<uint64_t, std::string>(
      kMetricComonentId, kMetricDesc, kMetricFieldOneDesc);
  auto* handle = metric->GetHandle(kMetricFieldStrValue);
  for (size_t i = 0; i < 1000000; ++i) {
    handle->AddValue(i);
  }
  metric_manager_.reset();

  NumericMetricsResultHandle* result_handle =
      MetricsParserParse(kTestOutput, ".*", kMetricFieldStrValue, 0, 0, 0);
  ASSERT_NE(nullptr, result_handle);
  ASSERT_FALSE(result_handle->Advance());

  result_handle = MetricsParserParse(kTestOutput, ".*", kMetricFieldStrValue,
                                     std::numeric_limits<uint64_t>::max(),
                                     std::numeric_limits<uint64_t>::max(), 0);
  ASSERT_NE(nullptr, result_handle);
  ASSERT_FALSE(result_handle->Advance());

  result_handle =
      MetricsParserParse(kTestOutput, ".*", kMetricFieldStrValue,
                         std::numeric_limits<uint64_t>::max(), 0, 0);
  ASSERT_NE(nullptr, result_handle);
  ASSERT_FALSE(result_handle->Advance());

  result_handle = MetricsParserParse(kTestOutput, ".*", kMetricFieldStrValue, 0,
                                     std::numeric_limits<uint64_t>::max(),
                                     std::numeric_limits<uint64_t>::max());
  ASSERT_NE(nullptr, result_handle);
  ASSERT_TRUE(result_handle->Advance());

  ASSERT_EQ(1ul, result_handle->Size());

  std::vector<double> values_out(1);
  std::vector<uint64_t> timestamps_out(1);
  MetricsParserResultHandleCopyInto(result_handle, &timestamps_out[0],
                                    &values_out[0]);
  ASSERT_EQ(999999, values_out.front());
  ASSERT_FALSE(result_handle->Advance());
  MetricsParserResultHandleFree(result_handle);
}

TEST_F(MetricFixture, ExternalMultiFieldsMatched) {
  auto* metric = metric_manager_->GetUnsafeMetric<uint64_t, std::string>(
      kMetricComonentId, kMetricDesc, kMetricFieldOneDesc);
  auto* handle_one = metric->GetHandle(kMetricFieldStrValue);
  auto* handle_two = metric->GetHandle(kMetricAnotherFieldStrValue);
  for (size_t i = 0; i < 1000000; ++i) {
    handle_one->AddValue(i + 1);
    handle_two->AddValue(i + 10);
  }
  metric_manager_.reset();

  NumericMetricsResultHandle* result_handle = MetricsParserParse(
      kTestOutput, ".*", ".*", 0, std::numeric_limits<uint64_t>::max(), 0);
  ASSERT_NE(nullptr, result_handle);

  std::map<std::pair<std::string, std::string>, std::vector<double>> values_map;
  while (result_handle->Advance()) {
    size_t size = MetricsParserResultHandleSize(result_handle);
    ASSERT_EQ(1000000ul, size);

    std::pair<std::string, std::string> key;
    key.first = result_handle->MetricString();
    key.second = result_handle->FieldString();

    std::vector<double> values_out(size);
    std::vector<uint64_t> timestamps_out(size);

    MetricsParserResultHandleCopyInto(result_handle, &timestamps_out[0],
                                      &values_out[0]);
    double sum = 0;
    for (double v : values_out) {
      sum += v;
    }
    values_map[key] = std::move(values_out);
  }
  MetricsParserResultHandleFree(result_handle);

  ASSERT_EQ(2ul, values_map.size());
  auto key_one = std::make_pair(kMetricComonentId, kMetricFieldStrValue);
  auto key_two = std::make_pair(kMetricComonentId, kMetricAnotherFieldStrValue);

  ASSERT_TRUE(ContainsKey(values_map, key_one));
  ASSERT_TRUE(ContainsKey(values_map, key_two));

  std::vector<double>& values_one = values_map[key_one];
  double total_one = std::accumulate(values_one.begin(), values_one.end(), 0.0);
  std::vector<double>& values_two = values_map[key_two];
  double total_two = std::accumulate(values_two.begin(), values_two.end(), 0.0);

  ASSERT_EQ((1000000.0 * 1000001.0) / 2.0, total_one);
  ASSERT_EQ((1000000.0 * 1000001.0) / 2.0 + 9000000.0, total_two);
}

}  // namespace test
}  // namespace metrics
}  // namespace ncode
