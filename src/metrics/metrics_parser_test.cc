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

#define CHECK_BAD_FIELDS(x) EXPECT_DEATH(FieldsMatcher::FromString(x), ".*");

#define CHECK_GOOD_FIELDS(x) \
  FieldsMatcher::FromString(Substitute(x, FieldsMatcher::kIntMatcher, 50));

#define CHECK_BAD_FIELDS_MULTI(x)                                             \
  EXPECT_DEATH(                                                               \
      FieldsMatcher::FromString(Substitute(x, FieldsMatcher::kIntMatcher, 50, \
                                           FieldsMatcher::kIntMatcher, 50)),  \
      ".*");

#define CHECK_GOOD_FIELDS_MULTI(x)                                        \
  FieldsMatcher::FromString(Substitute(x, FieldsMatcher::kIntMatcher, 50, \
                                       FieldsMatcher::kIntMatcher, 50));

// A bunch of bad fields strings.
TEST(FieldsMatcherString, BadStrings) {
  CHECK_BAD_FIELDS("");
  CHECK_BAD_FIELDS("junk");
  CHECK_BAD_FIELDS("junk(");
  CHECK_BAD_FIELDS("junk()");
  CHECK_BAD_FIELDS("junk( )");
  CHECK_BAD_FIELDS(Substitute("$0", FieldsMatcher::kStringMatcher));
  CHECK_BAD_FIELDS(Substitute("$0(", FieldsMatcher::kStringMatcher));
  CHECK_BAD_FIELDS(Substitute("$0)", FieldsMatcher::kStringMatcher));
  CHECK_BAD_FIELDS(Substitute("$0($1))", FieldsMatcher::kIntMatcher, 50));
  CHECK_BAD_FIELDS(Substitute("$0(($1)", FieldsMatcher::kIntMatcher, 50));
  CHECK_BAD_FIELDS(Substitute("($0($1)", FieldsMatcher::kIntMatcher, 50));
  CHECK_BAD_FIELDS(Substitute("($0($1))", FieldsMatcher::kIntMatcher, 50));
}

TEST(FieldsMatcherString, GoodStrings) {
  CHECK_GOOD_FIELDS("$0($1)");
  CHECK_GOOD_FIELDS("$0( $1)");
  CHECK_GOOD_FIELDS("$0( $1 )");
  CHECK_GOOD_FIELDS("$0( $1  )");
  CHECK_GOOD_FIELDS("$0 ($1)");
  CHECK_GOOD_FIELDS("$0 ( $1)");
  CHECK_GOOD_FIELDS("$0 ( $1 ) ");
}

TEST(FieldsMatcherString, BadStringsMulti) {
  CHECK_BAD_FIELDS_MULTI("$0($1)|$2($3)");
  CHECK_BAD_FIELDS_MULTI("$0($1),$2($3)");
}

TEST(FieldsMatcherString, GoodStringsMulti) {
  CHECK_GOOD_FIELDS_MULTI("$0($1)$2($3)");
  CHECK_GOOD_FIELDS_MULTI("$0($1) $2($3)");
  CHECK_GOOD_FIELDS_MULTI(" $0($1)   $2($3)  ");
}

static void CheckIntField(const FieldsMatcher& matcher, uint32_t value,
                          bool ok) {
  FieldsMatcher::FieldList field_list;
  PBMetricField* field = field_list.Add();
  field->set_type(PBMetricField::UINT32);
  field->set_uint32_value(value);
  if (ok) {
    ASSERT_TRUE(matcher.Matches(field_list));
  } else {
    ASSERT_FALSE(matcher.Matches(field_list));
  }

  field_list.Clear();
  field = field_list.Add();
  field->set_type(PBMetricField::UINT64);
  field->set_uint64_value(value);
  if (ok) {
    ASSERT_TRUE(matcher.Matches(field_list));
  } else {
    ASSERT_FALSE(matcher.Matches(field_list));
  }
}

static void CheckStringField(const FieldsMatcher& matcher,
                             const std::string& value, bool ok) {
  FieldsMatcher::FieldList field_list;
  PBMetricField* field = field_list.Add();
  field->set_type(PBMetricField::STRING);
  field->set_string_value(value);
  if (ok) {
    ASSERT_TRUE(matcher.Matches(field_list));
  } else {
    ASSERT_FALSE(matcher.Matches(field_list));
  }
}

TEST(FieldsMatcherString, SingleInt) {
  FieldsMatcher fields_matcher = FieldsMatcher::FromString(
      Substitute("$0($1)", FieldsMatcher::kIntMatcher, 50));

  FieldsMatcher::FieldList field_list;
  ASSERT_TRUE(fields_matcher.Matches(field_list));

  PBMetricField* field = field_list.Add();
  field->set_type(PBMetricField::BOOL);
  field->set_bool_value(true);
  ASSERT_FALSE(fields_matcher.Matches(field_list));

  CheckStringField(fields_matcher, "50", false);
  CheckIntField(fields_matcher, 49, false);
  CheckIntField(fields_matcher, 50, true);
  CheckIntField(fields_matcher, 500, false);

  field_list.Clear();
  field = field_list.Add();
  field->set_type(PBMetricField::UINT32);
  field->set_uint32_value(50);
  field = field_list.Add();
  field->set_type(PBMetricField::UINT64);
  field->set_uint64_value(500);
  ASSERT_TRUE(fields_matcher.Matches(field_list));

  field_list.Clear();
  field = field_list.Add();
  field->set_type(PBMetricField::UINT32);
  field->set_uint32_value(500);
  field = field_list.Add();
  field->set_type(PBMetricField::UINT64);
  field->set_uint64_value(50);
  ASSERT_FALSE(fields_matcher.Matches(field_list));
}

TEST(FieldsMatcherString, SingleIntLt) {
  FieldsMatcher fields_matcher = FieldsMatcher::FromString(
      Substitute("$0($1)", FieldsMatcher::kLtMatcher, 50));

  CheckIntField(fields_matcher, 50, false);
  CheckIntField(fields_matcher, 51, false);
  CheckIntField(fields_matcher, 49, true);
}

TEST(FieldsMatcherString, Combined) {
  FieldsMatcher matcher = FieldsMatcher::FromString(
      Substitute("$0($1) $2($3)", FieldsMatcher::kIntMatcher, 50,
                 FieldsMatcher::kStringMatcher, "some_string"));

  FieldsMatcher::FieldList field_list;
  PBMetricField* field = field_list.Add();
  field->set_type(PBMetricField::UINT32);
  field->set_uint32_value(50);

  field = field_list.Add();
  field->set_type(PBMetricField::STRING);
  field->set_string_value("some_string");

  ASSERT_TRUE(matcher.Matches(field_list));
}

TEST(FieldsMatcherString, SingleIntGt) {
  FieldsMatcher fields_matcher = FieldsMatcher::FromString(
      Substitute("$0($1)", FieldsMatcher::kGtMatcher, 50));

  CheckIntField(fields_matcher, 50, false);
  CheckIntField(fields_matcher, 51, true);
  CheckIntField(fields_matcher, 49, false);
}

TEST(FieldsMatcherString, SingleStringExact) {
  FieldsMatcher fields_matcher = FieldsMatcher::FromString(
      Substitute("$0($1)", FieldsMatcher::kStringMatcher, "some string"));

  CheckStringField(fields_matcher, "some string", true);
  CheckStringField(fields_matcher, "some other string", false);

  // The string is a regex, but it is more convenient if the beginning and end
  // of string characters are implied (anchored regex).
  CheckStringField(fields_matcher, "some string other", false);
  CheckStringField(fields_matcher, "other some string", false);
}

TEST(FieldsMatcherString, SingleStringRegex) {
  FieldsMatcher fields_matcher = FieldsMatcher::FromString(
      Substitute("$0($1)", FieldsMatcher::kStringMatcher, ".*"));

  CheckStringField(fields_matcher, "", true);
  CheckStringField(fields_matcher, "abracadabra", true);
}

TEST(FieldsMatcherString, SingleStringRegexNested) {
  FieldsMatcher fields_matcher = FieldsMatcher::FromString(Substitute(
      "$0($1)", FieldsMatcher::kStringMatcher, "abra(c(a)d)[ab]ra(\\.\\*)?"));

  CheckStringField(fields_matcher, "", false);
  CheckStringField(fields_matcher, "abracadara.*", true);
  CheckStringField(fields_matcher, "abracadara", true);
  CheckStringField(fields_matcher, "abracadabra", false);
}

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

  auto fields_matcher =
      FieldsMatcher::FromString(Substitute("string($0)", kMetricFieldStrValue));

  auto double_processor = make_unique<DoubleProcessor>(
      ".*", std::move(fields_matcher), double_callback);

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
      kTestOutput, ".*", Substitute("string($0)", kJunk).c_str(), 0,
      std::numeric_limits<uint64_t>::max(), 0);
  ASSERT_NE(nullptr, result_handle);
  ASSERT_FALSE(result_handle->Advance());
  MetricsParserResultHandleFree(result_handle);

  result_handle =
      MetricsParserParse(kTestOutput, kJunk,
                         Substitute("string($0)", kMetricFieldStrValue).c_str(),
                         0, std::numeric_limits<uint64_t>::max(), 0);
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
      MetricsParserParse(kTestOutput, kMetricComonentId,
                         Substitute("string($0)", kMetricFieldStrValue).c_str(),
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

  NumericMetricsResultHandle* numeric_result_handle = MetricsParserParse(
      kTestOutput, ".*", Substitute("string($0)", kMetricFieldStrValue).c_str(),
      0, std::numeric_limits<uint64_t>::max(), 0);
  ASSERT_FALSE(numeric_result_handle->Advance());

  BytesMetricsResultHandle* bytes_result_handle = MetricsParserBytesParse(
      kTestOutput, ".*", Substitute("string($0)", kMetricFieldStrValue).c_str(),
      0, std::numeric_limits<uint64_t>::max(), 0);
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

  NumericMetricsResultHandle* result_handle = MetricsParserParse(
      kTestOutput, ".*", Substitute("string($0)", kMetricFieldStrValue).c_str(),
      0, 0, 0);
  ASSERT_NE(nullptr, result_handle);
  ASSERT_FALSE(result_handle->Advance());

  result_handle = MetricsParserParse(
      kTestOutput, ".*", Substitute("string($0)", kMetricFieldStrValue).c_str(),
      std::numeric_limits<uint64_t>::max(),
      std::numeric_limits<uint64_t>::max(), 0);
  ASSERT_NE(nullptr, result_handle);
  ASSERT_FALSE(result_handle->Advance());

  result_handle = MetricsParserParse(
      kTestOutput, ".*", Substitute("string($0)", kMetricFieldStrValue).c_str(),
      std::numeric_limits<uint64_t>::max(), 0, 0);
  ASSERT_NE(nullptr, result_handle);
  ASSERT_FALSE(result_handle->Advance());

  result_handle = MetricsParserParse(
      kTestOutput, ".*", Substitute("string($0)", kMetricFieldStrValue).c_str(),
      0, std::numeric_limits<uint64_t>::max(),
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

  NumericMetricsResultHandle* result_handle =
      MetricsParserParse(kTestOutput, ".*", "string(.*)", 0,
                         std::numeric_limits<uint64_t>::max(), 0);
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
