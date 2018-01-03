#include "metrics.h"

#include <cstdint>
#include <limits>
#include <vector>
#include <thread>
#include <fcntl.h>
#include "metrics_test_util.h"
#include "metrics_parser.h"

#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include "gtest/gtest.h"
#include "ncode/substitute.h"

namespace nc {
namespace metrics {
namespace test {

TEST(Entry, Equality) {
  Entry<double> entry_one = {1.0, 5};
  Entry<double> entry_two = {1.0, 5};

  ASSERT_EQ(entry_one, entry_two);
  ASSERT_EQ(entry_two, entry_one);

  Entry<double> entry_three = {1.0, 4};
  ASSERT_NE(entry_one, entry_three);
  ASSERT_NE(entry_three, entry_one);

  Entry<double> entry_four = {2.0, 5};
  ASSERT_NE(entry_one, entry_four);
  ASSERT_NE(entry_four, entry_one);
}

// A metric with no fields.
TEST_F(MetricFixture, MetricNoFields) {
  auto* metric =
      metric_manager_->GetUnsafeMetric<uint64_t>(kMetricId, kMetricDesc);
  auto* handle_one = metric->GetHandle();
  auto* handle_two = metric->GetHandle();

  ASSERT_EQ(handle_one, handle_two);

  std::map<size_t, std::unique_ptr<PBManifestEntry>> index_to_manifest_entry =
      metric->ManifestIndexToManifestEntry();
  ASSERT_EQ(1ul, index_to_manifest_entry.size());

  const PBManifestEntry& manifest_entry = *index_to_manifest_entry[0];
  ASSERT_EQ(kMetricId, manifest_entry.id());
  ASSERT_EQ(kMetricDesc, manifest_entry.description());
  ASSERT_EQ(PBManifestEntry::UINT64, manifest_entry.type());
  ASSERT_EQ(0, manifest_entry.fields_size());
}

// This should create a single metric with one handle. The metric will store
// integer values and will have 2 fields -- one string and one integer.
TEST_F(MetricFixture, SingleHandle) {
  auto* metric =
      metric_manager_->GetUnsafeMetric<uint64_t, std::string, uint64_t>(
          kMetricId, kMetricDesc, kMetricFieldOneDesc, kMetricFieldTwoDesc);
  metric->GetHandle(kMetricFieldStrValue, kMetricFieldIntValue);

  std::map<size_t, std::unique_ptr<PBManifestEntry>> index_to_manifest_entry =
      metric->ManifestIndexToManifestEntry();
  ASSERT_EQ(1ul, index_to_manifest_entry.size());

  const PBManifestEntry& manifest_entry = *index_to_manifest_entry[0];
  ASSERT_EQ(kMetricId, manifest_entry.id());
  ASSERT_EQ(kMetricDesc, manifest_entry.description());
  ASSERT_EQ(PBManifestEntry::UINT64, manifest_entry.type());
  ASSERT_EQ(2, manifest_entry.fields_size());
  ASSERT_EQ(PBMetricField::STRING, manifest_entry.fields(0).type());
  ASSERT_EQ(kMetricFieldStrValue, manifest_entry.fields(0).string_value());
  ASSERT_EQ(PBMetricField::UINT64, manifest_entry.fields(1).type());
  ASSERT_EQ(kMetricFieldIntValue, manifest_entry.fields(1).uint64_value());
}

TEST_F(MetricFixture, SingleQuery) {
  auto* metric =
      metric_manager_->GetUnsafeMetric<double>(kMetricId, kMetricDesc);
  auto* handle = metric->GetHandle();
  handle->AddValue(5.0);

  Entry<double> entry;
  handle->MostRecentEntryInMem(&entry);

  ASSERT_EQ(5.0, entry.value);
  ASSERT_GE(timestamp_provider_->GetTimestamp(), entry.timestamp);
}

TEST_F(MetricFixture, MultiQuery) {
  auto* metric =
      metric_manager_->GetUnsafeMetric<double, std::string, uint64_t>(
          kMetricId, kMetricDesc, kMetricFieldOneDesc, kMetricFieldTwoDesc);
  auto* handle = metric->GetHandle(kMetricFieldStrValue, kMetricFieldIntValue);

  for (size_t i = 0; i < 10; ++i) {
    handle->AddValue(5.0 + i);
  }

  Entry<double> entry;
  ASSERT_EQ(true, handle->MostRecentEntryInMem(&entry));
  ASSERT_EQ(5.0 + 9, entry.value);
  ASSERT_GE(timestamp_provider_->GetTimestamp(), entry.timestamp);
}

static void ProcessEntriesFromDir(const std::string& metrics_dir,
                                  const std::string& metric_file,
                                  std::function<void(const PBMetricEntry&)> f) {
  using namespace google::protobuf;
  std::string file = nc::StrCat(metrics_dir, "/", metric_file);

  int fd = open(file.c_str(), O_RDONLY);
  CHECK(fd > 0);
  google::protobuf::io::FileInputStream input_stream(fd);

  uint32_t manifest_index;
  PBMetricEntry entry;
  while (true) {
    google::protobuf::io::CodedInputStream coded_stream(&input_stream);
    if (!metrics::parser::ReadDelimitedHeaderFrom(&manifest_index,
                                                  &coded_stream)) {
      break;
    }

    CHECK(metrics::parser::ReadDelimitedFrom(&entry, &coded_stream) == true);
    f(entry);
  }
}

TEST_F(MetricFixture, FileOutput) {
  auto* metric = metric_manager_->GetUnsafeMetric<double, std::string>(
      kMetricId, kMetricDesc, kMetricFieldOneDesc);
  auto* handle = metric->GetHandle(kMetricFieldStrValue);

  for (size_t i = 0; i < 1000000; ++i) {
    handle->AddValue(i);
  }

  metric_manager_.reset();

  // The file should contain 1M entries.
  std::vector<PBMetricEntry> all_entries;
  all_entries.reserve(1000000);
  ProcessEntriesFromDir(kTestOutput, kMetricId,
                        [&all_entries](const PBMetricEntry& entry) {
                          all_entries.emplace_back(entry);
                        });

  // 1M + 1 for the manifest entry.
  ASSERT_EQ(1000001ul, all_entries.size());

  // The first entry should be the manifest entry.
  const PBMetricEntry& first_entry = all_entries.front();
  ASSERT_TRUE(first_entry.has_manifest_entry());
  ASSERT_TRUE(first_entry.manifest_entry().type() == PBManifestEntry::DOUBLE);

  for (size_t i = 0; i < 1000000; ++i) {
    const PBMetricEntry& entry = all_entries[i + 1];
    ASSERT_EQ(i, entry.double_value());
  }
}

TEST_F(MetricFixture, FileOutputCallback) {
  auto* metric = metric_manager_->GetUnsafeMetric<double, std::string>(
      kMetricId, kMetricDesc, kMetricFieldOneDesc);

  size_t i = 0;
  metric->GetHandle([&i] {
    ++i;
    return i;
  }, kMetricFieldStrValue);

  for (size_t i = 0; i < 1000000; ++i) {
    metric->PollAllHandles();
  }

  metric_manager_.reset();

  // The file should contain 1M entries.
  std::vector<PBMetricEntry> all_entries;
  all_entries.reserve(1000000);
  ProcessEntriesFromDir(kTestOutput, kMetricId,
                        [&all_entries](const PBMetricEntry& entry) {
                          all_entries.emplace_back(entry);
                        });

  // 1M + 1 for the manifest entry.
  ASSERT_EQ(1000001ul, all_entries.size());
}

TEST_F(MetricFixture, BytesBlob) {
  std::string metric_file = std::string(kTestOutput);
  auto* metric = metric_manager_->GetUnsafeMetric<BytesBlob, std::string>(
      kMetricId, kMetricDesc, kMetricFieldOneDesc);
  auto* handle = metric->GetHandle(kMetricFieldStrValue);

  std::string value = "1234567890";

  BytesBlob blob;
  blob.set_bytes_value(value);
  handle->AddValue(blob);

  metric_manager_.reset();

  std::vector<PBMetricEntry> all_entries;
  ProcessEntriesFromDir(kTestOutput, kMetricId,
                        [&all_entries](const PBMetricEntry& entry) {
                          all_entries.emplace_back(entry);
                        });

  ASSERT_EQ(2ul, all_entries.size());
  const PBMetricEntry& entry = all_entries[1];
  ASSERT_EQ(value, entry.bytes_value().bytes_value());
}

TEST_F(MetricFixture, MultiThreadWriteToHandle) {
  auto* metric = metric_manager_->GetThreadSafeMetric<double, std::string>(
      kMetricId, kMetricDesc, kMetricFieldOneDesc);
  auto* handle = metric->GetHandle(kMetricFieldStrValue);

  // Writing to the same handle from different threads should preserve the
  // values.
  std::thread t1([handle] {
    for (size_t i = 0; i < 1000000; ++i) {
      handle->AddValue(5.0 + i);
    }
  });

  std::thread t2([handle] {
    for (size_t i = 0; i < 1000000; ++i) {
      handle->AddValue(5.0 + i);
    }
  });

  t1.join();
  t2.join();

  metric_manager_.reset();

  // The file should contain 2M entries, with 2 copies for each value.
  std::map<double, int> values;
  ProcessEntriesFromDir(kTestOutput, kMetricId,
                        [&values](const PBMetricEntry& entry) {
                          values[entry.double_value()] += 1;
                        });

  ASSERT_EQ(1000001ul, values.size());
  for (size_t i = 1; i < 1000001; ++i) {
    ASSERT_NE(values.end(), values.find(5.0 + i - 1));
    ASSERT_EQ(2, values[5.0 + i - 1]);
  }
}

TEST_F(MetricFixture, MultiThreadGetHandle) {
  std::string metric_file = std::string(kTestOutput);
  auto* metric = metric_manager_->GetThreadSafeMetric<double, std::string>(
      kMetricId, kMetricDesc, kMetricFieldOneDesc);
  std::array<std::vector<MetricHandle<double, true>*>, 2> ptr_storage;

  // Writing to the same handle from different threads should preserve the
  // values.
  std::thread t1([metric, &ptr_storage] {
    for (size_t i = 0; i < 1000000; ++i) {
      auto* handle =
          metric->GetHandle(Substitute("$0_$1", kMetricFieldStrValue, i));
      ptr_storage[0].emplace_back(handle);
    }
  });

  std::thread t2([metric, &ptr_storage] {
    for (size_t i = 0; i < 1000000; ++i) {
      auto* handle =
          metric->GetHandle(Substitute("$0_$1", kMetricFieldStrValue, i));
      ptr_storage[1].emplace_back(handle);
    }
  });

  t1.join();
  t2.join();

  ASSERT_EQ(1000000ul, ptr_storage[0].size());
  ASSERT_EQ(1000000ul, ptr_storage[1].size());
  for (size_t i = 0; i < 1000000ul; ++i) {
    ASSERT_EQ(ptr_storage[0][i], ptr_storage[1][i]);
  }
}

TEST_F(MetricFixture, PerMetricFileOutput) {
  std::string metric_id_one = "metric_one";
  std::string metric_id_two = "metric_two";

  auto* metric_one = metric_manager_->GetUnsafeMetric<double, std::string>(
      metric_id_one, kMetricDesc, kMetricFieldOneDesc);
  auto* metric_two = metric_manager_->GetUnsafeMetric<double, std::string>(
      metric_id_two, kMetricDesc, kMetricFieldOneDesc);

  auto* handle_one = metric_one->GetHandle(kMetricFieldStrValue);
  auto* handle_two = metric_two->GetHandle(kMetricFieldStrValue);

  handle_one->AddValue(10.0);
  handle_two->AddValue(20.0);

  metric_manager_.reset();

  std::vector<PBMetricEntry> entries_one;
  std::vector<PBMetricEntry> entries_two;
  ProcessEntriesFromDir(kTestOutput, metric_id_one,
                        [&entries_one](const PBMetricEntry& entry) {
                          entries_one.emplace_back(entry);
                        });
  ProcessEntriesFromDir(kTestOutput, metric_id_two,
                        [&entries_two](const PBMetricEntry& entry) {
                          entries_two.emplace_back(entry);
                        });

  ASSERT_EQ(2ul, entries_one.size());
  ASSERT_EQ(2ul, entries_two.size());
  ASSERT_EQ(20.0, entries_two.back().double_value());
  ASSERT_EQ(10.0, entries_one.back().double_value());
}

}  // namespace test
}  // namespace metrics
}  // namespace ncode
