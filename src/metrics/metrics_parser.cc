#include "metrics_parser.h"

#include <fcntl.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>

#include "ncode_common/src/strutil.h"
#include "ncode_common/src/substitute.h"

namespace nc {
namespace metrics {
namespace parser {

constexpr size_t kNumericFieldWidth = 10;
constexpr size_t kLongTextWidth = 100;
constexpr size_t kShortTextWidth = 40;

static std::string SingleFieldToString(const PBMetricField& field) {
  switch (field.type()) {
    case PBMetricField::BOOL:
      return std::to_string(field.bool_value());
    case PBMetricField::STRING:
      return field.string_value();
    case PBMetricField::UINT32:
      return std::to_string(field.uint32_value());
    case PBMetricField::UINT64:
      return std::to_string(field.uint64_value());
    default:
      LOG(FATAL) << "Don't know what to do with "
                 << PBMetricField_Type_Name(field.type());
  }
  return "";
}

std::string GetFieldString(const PBManifestEntry& entry) {
  std::string out;

  for (int i = 0; i < entry.fields_size(); ++i) {
    const PBMetricField& field = entry.fields(i);
    out += SingleFieldToString(field);
    if (i != entry.fields_size() - 1) {
      out += ':';
    }
  }

  if (out.empty()) {
    out = "EMPTY";
  }

  return out;
}

template <>
void ParseEntryFromProtobuf<uint64_t>(const PBMetricEntry& entry,
                                      Entry<uint64_t>* out) {
  out->timestamp = entry.timestamp();
  out->value = entry.uint64_value();
}

template <>
void ParseEntryFromProtobuf<uint32_t>(const PBMetricEntry& entry,
                                      Entry<uint32_t>* out) {
  out->timestamp = entry.timestamp();
  out->value = entry.uint32_value();
}

template <>
void ParseEntryFromProtobuf<bool>(const PBMetricEntry& entry,
                                  Entry<bool>* out) {
  out->timestamp = entry.timestamp();
  out->value = entry.bool_value();
}

template <>
void ParseEntryFromProtobuf<double>(const PBMetricEntry& entry,
                                    Entry<double>* out) {
  out->timestamp = entry.timestamp();
  out->value = entry.double_value();
}

template <>
void ParseEntryFromProtobuf<std::string>(const PBMetricEntry& entry,
                                         Entry<std::string>* out) {
  out->timestamp = entry.timestamp();
  out->value = entry.string_value();
}

template <>
void ParseEntryFromProtobuf<BytesBlob>(const PBMetricEntry& entry,
                                       Entry<BytesBlob>* out) {
  out->timestamp = entry.timestamp();
  out->value = entry.bytes_value();
}

InputStream::InputStream(const std::string& file) {
  fd_ = open(file.c_str(), O_RDONLY);
  CHECK(fd_ > 0) << "Bad input file " << file << ": " << strerror(errno);
  file_input_ = make_unique<google::protobuf::io::FileInputStream>(fd_);
  input_ =
      make_unique<google::protobuf::io::CodedInputStream>(file_input_.get());
  input_->SetTotalBytesLimit(std::numeric_limits<int>::max(),
                             std::numeric_limits<int>::max());
}

InputStream::~InputStream() {
  // close streams
  input_.reset();
  file_input_->Close();
  file_input_.reset();
  close(fd_);
}

bool InputStream::ReadDelimitedHeaderFrom(uint32_t* manifest_index) {
  // Read the manifest.
  if (!input_->ReadVarint32(manifest_index)) {
    return false;
  }
  return true;
}

bool InputStream::SkipMessage() {
  // Read the size.
  uint32_t size;
  if (!input_->ReadVarint32(&size)) {
    return false;
  }

  return input_->Skip(size);
}

bool InputStream::ReadDelimitedFrom(PBMetricEntry* message) {
  // Read the size.
  uint32_t size;
  if (!input_->ReadVarint32(&size)) {
    return false;
  }

  // Tell the stream not to read beyond that size.
  google::protobuf::io::CodedInputStream::Limit limit = input_->PushLimit(size);

  // Parse the message.
  if (!message->MergeFromCodedStream(input_.get())) {
    return false;
  }

  if (!input_->ConsumedEntireMessage()) {
    return false;
  }

  // Release the limit.
  input_->PopLimit(limit);

  return true;
}

void MetricsParser::Parse() {
  InputStream input_stream(metrics_file_);

  // This vector contains a map from the indices of manifests seen so far to the
  // list of processors that are interested in those manifests.
  std::vector<std::vector<MetricProcessor*>> manifest_index_to_processors;

  // Manifest entries.
  std::vector<std::unique_ptr<PBManifestEntry>> manifest_entries;

  uint32_t manifest_index;
  PBMetricEntry entry;
  while (true) {
    if (!input_stream.ReadDelimitedHeaderFrom(&manifest_index)) {
      break;
    }

    if (manifest_index == MetricBase::kManifestEntryMetaIndex) {
      // The following entry contains a manifest entry. Will read it in and pass
      // it to all processors to see if anyone is interested. A new vector will
      // be added to the back of manifest_index_to_processors and all interested
      // processors will be added to it.

      if (!input_stream.ReadDelimitedFrom(&entry)) {
        LOG(INFO) << "Unable to read in manifest entry";
        break;
      }

      CHECK(entry.has_manifest_entry())
          << "Wrong manifest index for manifest entry";

      std::vector<MetricProcessor*> interested;
      for (const auto& processor : processors_) {
        if (processor->InterestedIn(entry.manifest_entry(),
                                    manifest_index_to_processors.size())) {
          interested.emplace_back(processor.get());
        }
      }

      manifest_index_to_processors.emplace_back(interested);
      auto manifest_entry_ptr =
          std::unique_ptr<PBManifestEntry>(entry.release_manifest_entry());
      manifest_entries.emplace_back(std::move(manifest_entry_ptr));
      continue;
    }

    CHECK(manifest_index <= manifest_index_to_processors.size())
        << "Unknown manifest index";

    const std::vector<MetricProcessor*>& interested_processors =
        manifest_index_to_processors[manifest_index];
    if (interested_processors.empty()) {
      // No one is interested
      if (!input_stream.SkipMessage()) {
        LOG(INFO) << "Unable to skip entry";
        break;
      }
      continue;
    }

    // Have to read in the entire entry.
    if (!input_stream.ReadDelimitedFrom(&entry)) {
      LOG(INFO) << "Unable to read entry";
    }

    for (MetricProcessor* processor : interested_processors) {
      processor->ProcessEntry(entry, *manifest_entries[manifest_index],
                              manifest_index);
    }

    entry.Clear();
  }
}

static bool IsNumeric(const WrappedEntry& wrapped_entry) {
  PBManifestEntry::Type type = wrapped_entry.manifest_entry().type();
  return type == PBManifestEntry::DOUBLE || type == PBManifestEntry::UINT32 ||
         type == PBManifestEntry::UINT64;
}

static double ExtractNumericValueOrDie(PBManifestEntry::Type type,
                                       const PBMetricEntry& entry) {
  if (type == PBManifestEntry::DOUBLE) {
    return entry.double_value();
  }

  if (type == PBManifestEntry::UINT32) {
    return entry.uint32_value();
  }

  if (type == PBManifestEntry::UINT64) {
    return entry.uint64_value();
  }

  LOG(FATAL) << "No numeric value in entry";
  return 0;
}

void WrappedEntry::ChildEntry(bool numeric, double value) {
  ++num_entries_;
  if (numeric) {
    summary_stats_.Add(value);
    if (value > 0) {
      ++num_non_zero_entries_;
    }
  }
}

Manifest MetricsParser::ParseManifest() const {
  InputStream input_stream(metrics_file_);

  // Manifest entries.
  std::vector<std::unique_ptr<WrappedEntry>> all_entries;

  uint32_t manifest_index;
  PBMetricEntry entry;
  while (true) {
    if (!input_stream.ReadDelimitedHeaderFrom(&manifest_index)) {
      break;
    }

    if (manifest_index == MetricBase::kManifestEntryMetaIndex) {
      // The following entry contains a manifest entry.
      CHECK(input_stream.ReadDelimitedFrom(&entry))
          << "Unable to read in manifest entry";
      CHECK(entry.has_manifest_entry())
          << "Wrong manifest index for manifest entry";

      auto manifest_entry_ptr =
          std::unique_ptr<PBManifestEntry>(entry.release_manifest_entry());
      auto new_entry_ptr = nc::make_unique<WrappedEntry>(
          all_entries.size(), std::move(manifest_entry_ptr));
      all_entries.emplace_back(std::move(new_entry_ptr));
    } else {
      CHECK(manifest_index < all_entries.size())
          << "Unknown manifest index " << manifest_index
          << " only know indices up to " << all_entries.size();
      WrappedEntry& wrapped_entry = *(all_entries[manifest_index]);

      bool numeric;
      double value = 0.0;
      if ((numeric = IsNumeric(wrapped_entry))) {
        if (!input_stream.ReadDelimitedFrom(&entry)) {
          LOG(ERROR) << "Unable to read entry";
          break;
        }
        value = ExtractNumericValueOrDie(wrapped_entry.manifest_entry().type(),
                                         entry);
      } else {
        if (!input_stream.SkipMessage()) {
          LOG(ERROR) << "Unable to skip entry";
          break;
        }
      }

      wrapped_entry.ChildEntry(numeric, value);
    }
  }

  std::map<std::string, std::vector<std::unique_ptr<WrappedEntry>>>
      id_to_manifest;
  for (auto& wrapped_entry : all_entries) {
    const std::string& id = wrapped_entry->manifest_entry().id();
    id_to_manifest[id].emplace_back(std::move(wrapped_entry));
  }

  return {std::move(id_to_manifest)};
}

static constexpr char kMetricIdColumnName[] = "Metric Id";
static constexpr char kTypeColumnName[] = "Type";
static constexpr char kFieldsColumnName[] = "Fields";
static constexpr char kSetsCountColumnName[] = "Sets";
static constexpr char kValuesCountColumnName[] = "Values";
static constexpr char kNoFields[] = "NO FIELDS";

std::string Manifest::FullToString() const {
  std::stringstream ss;
  ss << std::setw(kShortTextWidth) << std::left << kMetricIdColumnName
     << std::setw(kShortTextWidth) << std::left << kTypeColumnName
     << std::setw(kLongTextWidth) << std::left << kFieldsColumnName
     << std::setw(kNumericFieldWidth) << std::left << kSetsCountColumnName
     << std::setw(kNumericFieldWidth) << std::left << kValuesCountColumnName
     << std::endl;
  for (const auto& id_and_manifest_entries : entries_) {
    const std::string& id = id_and_manifest_entries.first;
    const std::vector<std::unique_ptr<WrappedEntry>>& entries =
        id_and_manifest_entries.second;

    size_t total_values = 0;
    for (const auto& entry : entries) {
      total_values += entry->num_entries();
    }

    // Entries with the same id will have the same fields. The values of those
    // fields will be different.
    const WrappedEntry* first_entry = entries.front().get();
    ss << std::setw(kShortTextWidth) << std::left << id;
    ss << std::setw(kShortTextWidth) << std::left
       << PBManifestEntry_Type_Name(first_entry->manifest_entry().type());

    std::vector<std::string> fields_strings;
    for (const PBMetricField& field : first_entry->manifest_entry().fields()) {
      fields_strings.emplace_back(StrCat(PBMetricField::Type_Name(field.type()),
                                         "(", field.description(), ")"));
    }
    if (fields_strings.empty()) {
      fields_strings.emplace_back(kNoFields);
    }

    ss << std::setw(kLongTextWidth) << std::left << Join(fields_strings, ",");
    ss << std::setw(kNumericFieldWidth) << std::left
       << std::to_string(entries.size());
    ss << std::setw(kNumericFieldWidth) << std::left
       << std::to_string(total_values);
    ss << std::endl;
  }
  return ss.str();
}

static std::string GetFieldsString(
    const google::protobuf::RepeatedPtrField<PBMetricField>& fields) {
  std::vector<std::string> fields_strings;
  for (const PBMetricField& field : fields) {
    fields_strings.emplace_back(field.description());
  }
  if (fields_strings.empty()) {
    fields_strings.emplace_back(kNoFields);
  }

  return Join(fields_strings, ",");
}

std::string Manifest::ToString(const std::string& metric_id) const {
  const std::vector<std::unique_ptr<WrappedEntry>>& metric_entries =
      FindOrDie(entries_, metric_id);

  CHECK(!metric_entries.empty());
  const PBManifestEntry& first_entry = metric_entries.front()->manifest_entry();

  std::stringstream ss;
  ss << std::setw(kLongTextWidth) << std::left
     << GetFieldsString(first_entry.fields()) << std::endl;

  for (const auto& wrapped_entry : metric_entries) {
    size_t count = wrapped_entry->num_entries();
    const PBManifestEntry& entry = wrapped_entry->manifest_entry();
    ss << std::setw(kLongTextWidth) << std::left << GetFieldString(entry);
    ss << std::setw(kNumericFieldWidth) << std::left << std::to_string(count)
       << std::endl;
  }

  return ss.str();
}

uint64_t Manifest::TotalEntryCount() const {
  uint64_t total = 0;
  for (const auto& id_and_entries : entries_) {
    const std::vector<std::unique_ptr<WrappedEntry>>& entries_for_metric =
        id_and_entries.second;
    for (const auto& entry : entries_for_metric) {
      total += entry->num_entries();
    }
  }

  return total;
}

MetricsParser::MetricsParser(const std::string& metrics_file)
    : metrics_file_(metrics_file) {}

void NumericMetricsResultHandle::CopyInto(uint64_t* timestamps_out,
                                          double* values_out) {
  assert(next_it_ != id_to_values_.end());
  std::vector<ValuesAndManifest<double>::TimestampAndValue>& values =
      next_it_->second.values;
  for (size_t i = 0; i < values.size(); ++i) {
    timestamps_out[i] = values[i].first;
    values_out[i] = values[i].second;
  }
}

void BytesMetricsResultHandle::CopyInto(uint64_t* timestamps_out,
                                        char* values_out) {
  assert(next_it_ != id_to_values_.end());
  std::vector<ValuesAndManifest<std::string>::TimestampAndValue>& values =
      next_it_->second.values;

  size_t offset = 0;
  for (size_t i = 0; i < values.size(); ++i) {
    timestamps_out[i] = values[i].first;
    const std::string& string_to_copy = values[i].second;
    size_t string_size = string_to_copy.size();

    string_to_copy.copy(values_out + offset, string_size);
    offset += string_size;
  }
}

uint64_t BytesMetricsResultHandle::BufferSize(size_t i) {
  assert(next_it_ != id_to_values_.end());
  std::vector<ValuesAndManifest<std::string>::TimestampAndValue>& values =
      next_it_->second.values;
  assert(i < values.size());
  return values[i].second.size();
}

std::map<std::pair<std::string, std::string>,
         std::vector<std::pair<uint64_t, double>>>
SimpleParseNumericData(const std::string& metrics_file,
                       const std::string& metric_regix,
                       const std::string& fields_to_match,
                       uint64_t min_timestamp, uint64_t max_timestamp,
                       uint64_t limiting_timestamp) {
  std::map<std::pair<std::string, std::string>,
           std::vector<std::pair<uint64_t, double>>> out;

  auto result_handle =
      std::unique_ptr<NumericMetricsResultHandle>(MetricsParserParse(
          metrics_file.c_str(), metric_regix.c_str(), fields_to_match.c_str(),
          min_timestamp, max_timestamp, limiting_timestamp));
  if (!result_handle) {
    return out;
  }

  while (result_handle->Advance()) {
    std::string metric_id = result_handle->MetricString();
    std::string fields = result_handle->FieldString();
    std::vector<std::pair<uint64_t, double>>& vector = out[{metric_id, fields}];
    vector = std::move(result_handle->MutableValues());
  }

  return out;
}

std::map<std::pair<std::string, std::string>,
         std::vector<std::pair<uint64_t, double>>>
SimpleParseNumericData(const std::string& metrics_file,
                       const std::set<uint32_t> ids, uint64_t min_timestamp,
                       uint64_t max_timestamp, uint64_t limiting_timestamp) {
  using DoubleProcessor = IdCallbackProcessor<double, PBManifestEntry::DOUBLE>;
  using Uint32Processor =
      IdCallbackProcessor<uint32_t, PBManifestEntry::UINT32>;
  using Uint64Processor =
      IdCallbackProcessor<uint64_t, PBManifestEntry::UINT64>;

  auto handle = make_unique<NumericMetricsResultHandle>();
  DoubleProcessor::Callback double_callback = [&handle, min_timestamp,
                                               max_timestamp,
                                               limiting_timestamp](
      const Entry<double>& entry, const PBManifestEntry& manifest_entry,
      uint32_t manifest_index) {
    if (entry.timestamp < max_timestamp && entry.timestamp >= min_timestamp) {
      handle->Update(entry.timestamp, entry.value, manifest_index,
                     manifest_entry, limiting_timestamp);
    }
  };

  Uint32Processor::Callback uint32_callback = [&handle, min_timestamp,
                                               max_timestamp,
                                               limiting_timestamp](
      const Entry<uint32_t>& entry, const PBManifestEntry& manifest_entry,
      uint32_t manifest_index) {
    if (entry.timestamp < max_timestamp && entry.timestamp >= min_timestamp) {
      handle->Update(entry.timestamp, static_cast<double>(entry.value),
                     manifest_index, manifest_entry, limiting_timestamp);
    }
  };

  Uint64Processor::Callback uint64_callback = [&handle, min_timestamp,
                                               max_timestamp,
                                               limiting_timestamp](
      const Entry<uint64_t>& entry, const PBManifestEntry& manifest_entry,
      int32_t manifest_index) {
    if (entry.timestamp < max_timestamp && entry.timestamp >= min_timestamp) {
      handle->Update(entry.timestamp, static_cast<double>(entry.value),
                     manifest_index, manifest_entry, limiting_timestamp);
    }
  };

  auto double_processor = make_unique<DoubleProcessor>(ids, double_callback);
  auto uint32_processor = make_unique<Uint32Processor>(ids, uint32_callback);
  auto uint64_processor = make_unique<Uint64Processor>(ids, uint64_callback);

  MetricsParser parser(metrics_file);
  parser.AddProcessor(std::move(double_processor));
  parser.AddProcessor(std::move(uint32_processor));
  parser.AddProcessor(std::move(uint64_processor));

  parser.Parse();
  handle->Sort();

  std::map<std::pair<std::string, std::string>,
           std::vector<std::pair<uint64_t, double>>> out;
  while (handle->Advance()) {
    std::string metric_id = handle->MetricString();
    std::string fields = handle->FieldString();
    std::vector<std::pair<uint64_t, double>>& vector = out[{metric_id, fields}];
    vector = std::move(handle->MutableValues());
  }

  return out;
}

NumericMetricsResultHandle* MetricsParserParse(const char* metrics_file,
                                               const char* metric_regex,
                                               const char* fields_to_match,
                                               uint64_t min_timestamp,
                                               uint64_t max_timestamp,
                                               uint64_t limiting_timestamp) {
  using DoubleProcessor =
      QueryCallbackProcessor<double, PBManifestEntry::DOUBLE>;
  using Uint32Processor =
      QueryCallbackProcessor<uint32_t, PBManifestEntry::UINT32>;
  using Uint64Processor =
      QueryCallbackProcessor<uint64_t, PBManifestEntry::UINT64>;

  NumericMetricsResultHandle* return_handle = new NumericMetricsResultHandle;

  DoubleProcessor::Callback double_callback = [return_handle, min_timestamp,
                                               max_timestamp,
                                               limiting_timestamp](
      const Entry<double>& entry, const PBManifestEntry& manifest_entry,
      uint32_t manifest_index) {
    if (entry.timestamp < max_timestamp && entry.timestamp >= min_timestamp) {
      return_handle->Update(entry.timestamp, entry.value, manifest_index,
                            manifest_entry, limiting_timestamp);
    }
  };

  Uint32Processor::Callback uint32_callback = [return_handle, min_timestamp,
                                               max_timestamp,
                                               limiting_timestamp](
      const Entry<uint32_t>& entry, const PBManifestEntry& manifest_entry,
      uint32_t manifest_index) {
    if (entry.timestamp < max_timestamp && entry.timestamp >= min_timestamp) {
      return_handle->Update(entry.timestamp, static_cast<double>(entry.value),
                            manifest_index, manifest_entry, limiting_timestamp);
    }
  };

  Uint64Processor::Callback uint64_callback = [return_handle, min_timestamp,
                                               max_timestamp,
                                               limiting_timestamp](
      const Entry<uint64_t>& entry, const PBManifestEntry& manifest_entry,
      int32_t manifest_index) {
    if (entry.timestamp < max_timestamp && entry.timestamp >= min_timestamp) {
      return_handle->Update(entry.timestamp, static_cast<double>(entry.value),
                            manifest_index, manifest_entry, limiting_timestamp);
    }
  };

  auto double_processor = make_unique<DoubleProcessor>(
      metric_regex, fields_to_match, double_callback);
  auto uint32_processor = make_unique<Uint32Processor>(
      metric_regex, fields_to_match, uint32_callback);
  auto uint64_processor = make_unique<Uint64Processor>(
      metric_regex, fields_to_match, uint64_callback);

  MetricsParser parser(metrics_file);
  parser.AddProcessor(std::move(double_processor));
  parser.AddProcessor(std::move(uint32_processor));
  parser.AddProcessor(std::move(uint64_processor));

  parser.Parse();
  return_handle->Sort();
  return return_handle;
}

BytesMetricsResultHandle* MetricsParserBytesParse(const char* metrics_file,
                                                  const char* metric_regex,
                                                  const char* fields_to_match,
                                                  uint64_t min_timestamp,
                                                  uint64_t max_timestamp,
                                                  uint64_t limiting_timestamp) {
  using BytesProcessor =
      QueryCallbackProcessor<BytesBlob, PBManifestEntry::BYTES>;

  BytesMetricsResultHandle* return_handle = new BytesMetricsResultHandle;
  BytesProcessor::Callback bytes_callback = [return_handle, min_timestamp,
                                             max_timestamp, limiting_timestamp](
      const Entry<BytesBlob>& entry, const PBManifestEntry& manifest_entry,
      uint32_t manifest_index) {
    if (entry.timestamp < max_timestamp && entry.timestamp >= min_timestamp) {
      return_handle->Update(entry.timestamp, entry.value.bytes_value(),
                            manifest_index, manifest_entry, limiting_timestamp);
    }
  };

  auto bytes_processor = make_unique<BytesProcessor>(
      metric_regex, fields_to_match, bytes_callback);

  MetricsParser parser(metrics_file);
  parser.AddProcessor(std::move(bytes_processor));

  parser.Parse();
  return_handle->Sort();
  return return_handle;
}

bool MetricsParserResultHandleAdvance(NumericMetricsResultHandle* handle) {
  return handle->Advance();
}

uint64_t MetricsParserResultHandleSize(NumericMetricsResultHandle* handle) {
  return handle->Size();
}

bool MetricsParserBytesResultHandleAdvance(BytesMetricsResultHandle* handle) {
  return handle->Advance();
}

uint64_t MetricsParserBytesResultHandleSize(BytesMetricsResultHandle* handle) {
  return handle->Size();
}

uint64_t MetricsParserBytesResultHandleBufferSize(
    BytesMetricsResultHandle* handle, uint64_t i) {
  return handle->BufferSize(i);
}

// Returns a newly malloc-ed c-like array copy of the given string.
static char* StringToCString(const std::string& str) {
  size_t size = (str.size() + 1) * sizeof(char);
  char* c_str = static_cast<char*>(malloc(size));
  std::strncpy(c_str, str.c_str(), str.size());
  c_str[str.size()] = '\0';
  return c_str;
}

char* MetricsParserManifestSummary(const char* metrics_file) {
  MetricsParser parser(metrics_file);
  return StringToCString(parser.ParseManifest().FullToString());
}

char* MetricsParserManifestMetricSummary(const char* metrics_file,
                                         const char* metric_id) {
  MetricsParser parser(metrics_file);
  return StringToCString(parser.ParseManifest().ToString(metric_id));
}

char* MetricsParserResultHandleFieldString(NumericMetricsResultHandle* handle) {
  std::string return_string = handle->FieldString();
  return StringToCString(return_string);
}

char* MetricsParserResultHandleMetricString(
    NumericMetricsResultHandle* handle) {
  std::string return_string = handle->MetricString();
  return StringToCString(return_string);
}

char* MetricsParserBytesResultHandleFieldString(
    NumericMetricsResultHandle* handle) {
  std::string return_string = handle->FieldString();
  return StringToCString(return_string);
}

char* MetricsParserBytesResultHandleMetricString(
    NumericMetricsResultHandle* handle) {
  std::string return_string = handle->MetricString();
  return StringToCString(return_string);
}

void MetricsParserResultHandleCopyInto(NumericMetricsResultHandle* handle,
                                       uint64_t* timestamps_out,
                                       double* values_out) {
  handle->CopyInto(timestamps_out, values_out);
}

void MetricsParserBytesResultHandleCopyInto(BytesMetricsResultHandle* handle,
                                            uint64_t* timestamps_out,
                                            char* values_out) {
  handle->CopyInto(timestamps_out, values_out);
}

void MetricsParserResultHandleFree(NumericMetricsResultHandle* handle) {
  delete handle;
}

void MetricsParserBytesResultHandleFree(BytesMetricsResultHandle* handle) {
  delete handle;
}

void MetricsParserStringFree(char* str) { delete str; }
}  // namespace parser
}  // namespace metrics
}  // namespace ncode
