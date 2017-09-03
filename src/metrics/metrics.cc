#include "metrics.h"

#include <google/protobuf/repeated_field.h>
#include <google/protobuf/io/coded_stream.h>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "gflags/gflags.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/strutil.h"

DEFINE_string(metrics_output, "metrics.out",
              "Directory where the metrics will be stored");

namespace nc {
namespace metrics {

template <>
void SaveEntryToProtobuf<uint64_t>(const Entry<uint64_t>& entry,
                                   PBMetricEntry* out) {
  if (entry.timestamp) {
    out->set_timestamp(entry.timestamp);
  }
  out->set_uint64_value(entry.value);
}

template <>
void SaveEntryToProtobuf<uint32_t>(const Entry<uint32_t>& entry,
                                   PBMetricEntry* out) {
  if (entry.timestamp) {
    out->set_timestamp(entry.timestamp);
  }
  out->set_uint32_value(entry.value);
}

template <>
void SaveEntryToProtobuf<bool>(const Entry<bool>& entry, PBMetricEntry* out) {
  if (entry.timestamp) {
    out->set_timestamp(entry.timestamp);
  }
  out->set_bool_value(entry.value);
}

template <>
void SaveEntryToProtobuf<double>(const Entry<double>& entry,
                                 PBMetricEntry* out) {
  if (entry.timestamp) {
    out->set_timestamp(entry.timestamp);
  }
  out->set_double_value(entry.value);
}

template <>
void SaveEntryToProtobuf<std::string>(const Entry<std::string>& entry,
                                      PBMetricEntry* out) {
  if (entry.timestamp) {
    out->set_timestamp(entry.timestamp);
  }
  out->set_string_value(entry.value);
}

template <>
void SaveEntryToProtobuf<BytesBlob>(const Entry<BytesBlob>& entry,
                                    PBMetricEntry* out) {
  if (entry.timestamp) {
    out->set_timestamp(entry.timestamp);
  }
  *out->mutable_bytes_value() = entry.value;
}

OutputStream::OutputStream(const std::string& file) {
  fd_ = open(file.c_str(), O_WRONLY | O_TRUNC | O_CREAT,  // open mode
             S_IREAD | S_IWRITE | S_IRGRP | S_IROTH | S_ISUID);
  CHECK(fd_ > 0) << "Bad output file " << file;
  file_output_ = make_unique<google::protobuf::io::FileOutputStream>(fd_);
}

OutputStream::~OutputStream() {
  // close streams
  file_output_->Close();
  file_output_.reset();
  close(fd_);
}

void OutputStream::WriteBulk(const std::vector<PBMetricEntry>& entries,
                             uint32_t manifest_index) {
  std::lock_guard<std::mutex> lock(mu_);
  for (const auto& entry : entries) {
    CHECK(WriteDelimitedTo(entry, manifest_index));
  }
}

void OutputStream::WriteSingle(const PBMetricEntry& entry,
                               uint32_t manifest_index) {
  std::lock_guard<std::mutex> lock(mu_);
  CHECK(WriteDelimitedTo(entry, manifest_index));
}

// Writes a protobuf to the stream.
bool OutputStream::WriteDelimitedTo(const PBMetricEntry& entry,
                                    uint32_t manifest_index) {
  // Write the size and the manifest index. The index comes first.
  ::google::protobuf::io::CodedOutputStream coded_output(file_output_.get());

  coded_output.WriteVarint32(manifest_index);
  const int size = entry.ByteSize();
  coded_output.WriteVarint32(size);

  uint8_t* buffer = coded_output.GetDirectBufferForNBytesAndAdvance(size);
  if (buffer != nullptr) {
    // Optimization:  The message fits in one buffer, so use the faster
    // direct-to-array serialization path.
    entry.SerializeWithCachedSizesToArray(buffer);
  } else {
    // Slightly-slower path when the message is multiple buffers.
    entry.SerializeWithCachedSizes(&coded_output);
    if (coded_output.HadError()) {
      return false;
    }
  }
  return true;
}

void PopulateManifestEntryField(PBMetricField* field, uint64_t value) {
  field->set_type(PBMetricField::UINT64);
  field->set_uint64_value(value);
}

void PopulateManifestEntryField(PBMetricField* field, uint32_t value) {
  field->set_type(PBMetricField::UINT32);
  field->set_uint32_value(value);
}

void PopulateManifestEntryField(PBMetricField* field, bool value) {
  field->set_type(PBMetricField::BOOL);
  field->set_bool_value(value);
}

void PopulateManifestEntryField(PBMetricField* field,
                                const std::string& value) {
  field->set_type(PBMetricField::STRING);
  field->set_string_value(value);
}

MetricManager::MetricManager()
    : timestamp_provider_(make_unique<DefaultTimestampProvider>()) {}

void MetricManager::SetOutput(const std::string& output) {
  std::lock_guard<std::mutex> lock(mu_);
  if (!output_directory_.empty()) {
    LOG(ERROR) << "Output already set. Will ignore";
    return;
  }

  // If any of the metrics are already locked setting the stream will result in
  // broken output.
  for (const auto& metric : all_metrics_) {
    CHECK(!metric->stream_locked());
  }

  File::RecursivelyCreateDir(output, 0700);
  output_directory_ = output;
  for (const auto& metric : all_metrics_) {
    std::string local_output = StrCat(output_directory_, "/", metric->id());
    auto local_output_stream = make_unique<OutputStream>(local_output);
    metric->SetLocalOutputStream(std::move(local_output_stream));
  }
}

const TimestampProviderInterface* MetricManager::timestamp_provider() const {
  return timestamp_provider_.get();
}

const TimestampProviderInterface* MetricBase::timestamp_provider() const {
  return parent_manager_->timestamp_provider();
}

void MetricBase::SetLocalOutputStream(
    std::unique_ptr<OutputStream> output_stream) {
  CHECK(!local_output_stream_);
  CHECK(local_current_index_ == std::numeric_limits<size_t>::max());
  CHECK(!stream_locked_);
  local_output_stream_ = std::move(output_stream);
}

OutputStream* MetricBase::OutputStreamOrNull() {
  // Regardless of whether there was a stream or not we cannot allow further
  // changes to the stream, as it may result in missing manifest entries when
  // parsing.
  stream_locked_ = true;

  if (local_output_stream_) {
    return local_output_stream_.get();
  }

  return nullptr;
}

size_t MetricBase::NextIndex() { return ++local_current_index_; }

MetricHandleBase::MetricHandleBase(bool has_callback, MetricBase* parent_metric)
    : has_callback_(has_callback),
      metric_index_(0),
      parent_metric_(parent_metric) {}

std::string MetricHandleBase::TimestampToString(uint64_t timestamp) const {
  return parent_metric_->timestamp_provider()->TimestampToString(timestamp);
}

void MetricManager::PersistAllMetrics() {
  for (const auto& metric_ptr : all_metrics_) {
    metric_ptr->PersistAllHandles();
  }
}

void MetricManager::PollAllMetrics() {
  for (const auto& metric_ptr : all_metrics_) {
    metric_ptr->PollAllHandles();
  }
}

MetricManager::~MetricManager() { PersistAllMetrics(); }

MetricManager* DefaultMetricManager() {
  static MetricManager default_manager;
  return &default_manager;
}

void InitMetrics() {
  MetricManager* manager = DefaultMetricManager();
  manager->SetOutput(FLAGS_metrics_output);
}

void PopulateManifestEntryType(PBManifestEntry* out, uint64_t* dummy) {
  Unused(dummy);
  out->set_type(PBManifestEntry::UINT64);
}

void PopulateManifestEntryType(PBManifestEntry* out, uint32_t* dummy) {
  Unused(dummy);
  out->set_type(PBManifestEntry::UINT32);
}

void PopulateManifestEntryType(PBManifestEntry* out, bool* dummy) {
  Unused(dummy);
  out->set_type(PBManifestEntry::BOOL);
}

void PopulateManifestEntryType(PBManifestEntry* out, std::string* dummy) {
  Unused(dummy);
  out->set_type(PBManifestEntry::STRING);
}

void PopulateManifestEntryType(PBManifestEntry* out, double* dummy) {
  Unused(dummy);
  out->set_type(PBManifestEntry::DOUBLE);
}

void PopulateManifestEntryType(PBManifestEntry* out, BytesBlob* dummy) {
  Unused(dummy);
  out->set_type(PBManifestEntry::BYTES);
}

DefaultMetricManagerPoller::DefaultMetricManagerPoller(
    std::chrono::milliseconds period, EventQueue* event_queue)
    : EventConsumer("metric_poller", event_queue), period_(period) {
  EnqueueNext();
}

void DefaultMetricManagerPoller::HandleEvent() {
  metrics::DefaultMetricManager()->PollAllMetrics();
  EnqueueNext();
}

}  // namespace metrics
}  // namespace ncode
