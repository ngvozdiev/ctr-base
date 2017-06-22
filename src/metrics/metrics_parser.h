#ifndef METRICS_PARSER_H
#define METRICS_PARSER_H

#include <google/protobuf/repeated_field.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <cstdbool>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "metrics.pb.h"
#include "metrics.h"

namespace nc {
namespace metrics {
namespace parser {

template <typename T>
void ParseEntryFromProtobuf(const PBMetricEntry& entry, Entry<T>* out) {
  Unused(entry);
  Unused(out);
  assert(false);
}

template <>
void ParseEntryFromProtobuf<uint64_t>(const PBMetricEntry& entry,
                                      Entry<uint64_t>* out);

template <>
void ParseEntryFromProtobuf<uint32_t>(const PBMetricEntry& entry,
                                      Entry<uint32_t>* out);

template <>
void ParseEntryFromProtobuf<bool>(const PBMetricEntry& entry, Entry<bool>* out);

template <>
void ParseEntryFromProtobuf<double>(const PBMetricEntry& entry,
                                    Entry<double>* out);

template <>
void ParseEntryFromProtobuf<std::string>(const PBMetricEntry& entry,
                                         Entry<std::string>* out);

template <>
void ParseEntryFromProtobuf<BytesBlob>(const PBMetricEntry& entry,
                                       Entry<BytesBlob>* out);

class InputStream {
 public:
  InputStream(const std::string& file);

  ~InputStream();

  // Reads the header from a stream. This is usually followed by
  // ReadDelimitedFrom if the client thinks the manifest index is interesting or
  // calling SkipMessage and calling again if it's not.
  bool ReadDelimitedHeaderFrom(uint32_t* manifest_index);

  // If called after ReadDelimitedHeaderFrom will skip to the next message.
  bool SkipMessage();

  // Reads a protobuf from a stream. Should be called after
  // ReadDelimitedHeaderFrom.
  bool ReadDelimitedFrom(PBMetricEntry* message);

 private:
  // File descriptor. Closed on destruction.
  int fd_;

  // File input stream. Owned by this object.
  std::unique_ptr<google::protobuf::io::FileInputStream> file_input_;

  // Coded input stream.
  std::unique_ptr<google::protobuf::io::CodedInputStream> input_;

  DISALLOW_COPY_AND_ASSIGN(InputStream);
};

// Returns a human-readable string that described the fields in the entry.
std::string GetFieldString(const PBManifestEntry& entry);

// A class that knows how to match against a single field.
class SingleFieldMatcher {
 public:
  static std::unique_ptr<SingleFieldMatcher> FromString(
      const std::string& matcher_name, const std::string& matcher_string);

  virtual ~SingleFieldMatcher() {}

  virtual bool Matches(const PBMetricField& field) const = 0;
};

// Matches numeric values exactly.
class NumericValueExactMatcher : public SingleFieldMatcher {
 public:
  NumericValueExactMatcher(uint64_t value) : value_(value) {}

  bool Matches(const PBMetricField& field) const;

 private:
  uint64_t value_;
};

// Matches if a value is in range.
class NumericValueRangeMatcher : public SingleFieldMatcher {
 public:
  NumericValueRangeMatcher(uint64_t min, uint64_t max) : min_(min), max_(max) {}

  bool Matches(const PBMetricField& field) const;

 private:
  uint64_t min_;
  uint64_t max_;
};

// Matches string fields based on a regex. Case is ignored.
class StringRegexMatcher : public SingleFieldMatcher {
 public:
  StringRegexMatcher(const std::string& regex_string);

  bool Matches(const PBMetricField& field) const;

 private:
  std::regex field_regex_;
};

// A class that knows how to match against a set of fields.
class FieldsMatcher {
 public:
  using FieldList = ::google::protobuf::RepeatedPtrField<PBMetricField>;

  static constexpr char kIntMatcher[] = "int";
  static constexpr char kLtMatcher[] = "lt";
  static constexpr char kGtMatcher[] = "gt";
  static constexpr char kStringMatcher[] = "string";

  // Constructs a field matcher from a string. The string contains
  // comma-separated matcher descriptions. The format is as follows:
  // int(5) -- matches the integer value (5 in this case) exactly
  // int_lt(4) -- matches anything less than 4
  // int_gt(4) -- matches anything greater than 4
  // string(value) -- the field is a string that matches the value regex
  //
  // For example: int(5),string(ab[cd]) will match fields 5,abc and 5,abd
  static FieldsMatcher FromString(const std::string& str);

  // Similar to above, but instead of CHECK-failing on illegal syntax will
  // return false and not populate the output argument.
  static bool FromString(const std::string& str, FieldsMatcher* out);

  FieldsMatcher(FieldsMatcher&& other)
      : matchers_(std::move(other.matchers_)) {}

  FieldsMatcher(std::vector<std::unique_ptr<SingleFieldMatcher>> matchers)
      : matchers_(std::move(matchers)) {}

  // Returns true if the fields match.
  bool Matches(const FieldList& fields) const;

 private:
  // The matchers. All of them have to match incoming fields. The matcher at
  // index i has to match an incoming field at index i. If there are more
  // matchers than fields the extra matchers are ignored. If there are less
  // matchers than fields all extra fields are considered matched.
  std::vector<std::unique_ptr<SingleFieldMatcher>> matchers_;
};

// A generic interface for classes that know how to process metrics.
class MetricProcessor {
 public:
  // A processor that will express interest in metrics that match the .
  MetricProcessor() {}

  virtual ~MetricProcessor() {}

  // Returns true if this processor is interested in the metric identified by
  // the manifest entry.
  virtual bool InterestedIn(const PBManifestEntry& manifest_entry,
                            uint32_t manifest_index) = 0;

  // Processes an entry. This entry will come from one of the metrics that the
  // processor expressed interest in.
  virtual void ProcessEntry(const PBMetricEntry& entry,
                            const PBManifestEntry& manifest_entry,
                            uint32_t manifest_index) = 0;

 protected:
  DISALLOW_COPY_AND_ASSIGN(MetricProcessor);
};

// A class that processes a single metric entry and calls a callback with the
// timestamp and the value.
template <typename T>
class CallbackProcessor : public MetricProcessor {
 public:
  // A callback to call when an entry is processed. The first argument is the
  // entry, the second is the id, the third is the field string.
  typedef std::function<void(const Entry<T>&,
                             const PBManifestEntry& manifest_entry,
                             uint32_t manifest_index)> Callback;

  CallbackProcessor(Callback callback) : callback_(callback) {}

  void ProcessEntry(const PBMetricEntry& entry,
                    const PBManifestEntry& manifest_entry,
                    uint32_t manifest_index) override {
    ParseEntryFromProtobuf(entry, &scratch_entry_);
    callback_(scratch_entry_, manifest_entry, manifest_index);
  }

 private:
  // A callback to call when an entry is processed. The first argument is the
  // entry, the second is the id.
  Callback callback_;

  // A scratch entry to deserialize values into.
  Entry<T> scratch_entry_;
};

template <typename T, PBManifestEntry::Type EntryType>
class QueryCallbackProcessor : public CallbackProcessor<T> {
 public:
  QueryCallbackProcessor(const std::string& metric_regex_string,
                         FieldsMatcher fields_matcher,
                         typename CallbackProcessor<T>::Callback callback)
      : CallbackProcessor<T>(callback),
        original_string_(metric_regex_string),
        metric_regex_(metric_regex_string, std::regex_constants::icase),
        fields_matcher_(std::move(fields_matcher)) {}

  bool InterestedIn(const PBManifestEntry& manifest_entry,
                    uint32_t manifest_index) override {
    Unused(manifest_index);
    if (manifest_entry.type() != EntryType) {
      return false;
    }

    if (!std::regex_search(manifest_entry.id(), metric_regex_)) {
      return false;
    }

    return fields_matcher_.Matches(manifest_entry.fields());
  }

 private:
  std::string original_string_;

  // Regex the metric should match against.
  std::regex metric_regex_;

  // The fields of the metric will be matched against this.
  FieldsMatcher fields_matcher_;
};

// A callback processor that can be used if the manifest ids are known ahead of
// time.
template <typename T, PBManifestEntry::Type EntryType>
class IdCallbackProcessor : public CallbackProcessor<T> {
 public:
  IdCallbackProcessor(const std::set<uint32_t>& manifest_indices,
                      typename CallbackProcessor<T>::Callback callback)
      : CallbackProcessor<T>(callback), manifest_indices_(manifest_indices) {}

  bool InterestedIn(const PBManifestEntry& manifest_entry,
                    uint32_t manifest_index) override {
    if (manifest_entry.type() != EntryType) {
      return false;
    }
    return ContainsKey(manifest_indices_, manifest_index);
  }

 private:
  std::set<uint32_t> manifest_indices_;
};

// Combination of a manifest entry and the number of entries that belong to it
// in the entire file.
class WrappedEntry {
 public:
  WrappedEntry(size_t manifest_index,
               std::unique_ptr<PBManifestEntry> manifest_entry)
      : manifest_index_(manifest_index),
        manifest_entry_(std::move(manifest_entry)),
        num_entries_(0),
        num_non_zero_entries_(0),
        sum_(0) {}

  WrappedEntry(WrappedEntry&& other)
      : manifest_index_(other.manifest_index_),
        manifest_entry_(std::move(other.manifest_entry_)),
        num_entries_(other.num_entries_),
        num_non_zero_entries_(other.num_non_zero_entries_),
        sum_(other.sum_) {
    other.num_entries_ = 0;
    other.num_non_zero_entries_ = 0;
    other.sum_ = 0;
  }

  uint64_t manifest_index() const { return manifest_index_; }

  const PBManifestEntry& manifest_entry() const { return *manifest_entry_; }

  uint64_t num_entries() const { return num_entries_; }

  uint64_t num_non_zero_entries() const { return num_non_zero_entries_; }

  double sum() const { return sum_; }

  // Called for each entry that belongs to this manifest entry.
  void ChildEntry(bool numeric, double value);

 private:
  uint64_t manifest_index_;

  std::unique_ptr<PBManifestEntry> manifest_entry_;

  // Total number of entries.
  uint64_t num_entries_;

  // Number of non-zero entries. Only meaningful if the metric is numeric.
  uint64_t num_non_zero_entries_;

  // Total sum of entries. Only meaningful if the metric is numeric.
  double sum_;
};

class Manifest {
 public:
  Manifest(std::map<std::string, std::vector<WrappedEntry>>&& entries)
      : entries_(std::move(entries)) {}

  Manifest(Manifest&& other) : entries_(std::move(other.entries_)) {}

  // Returns a string summary of all metrics in the manifest.
  std::string FullToString() const;

  // Returns a string summary of the fields of the metric whose id exactly
  // matches the argument.
  std::string ToString(const std::string& metric_id) const;

  // The total number of entries across all child entries.
  uint64_t TotalEntryCount() const;

 private:
  // Entries, grouped by metric id.
  std::map<std::string, std::vector<WrappedEntry>> entries_;

  DISALLOW_COPY_AND_ASSIGN(Manifest);
};

// A class that can parse entries stored in a protobuf file.
class MetricsParser {
 public:
  MetricsParser(const std::string& metrics_file);

  // Adds a new processor to the parser.
  void AddProcessor(std::unique_ptr<MetricProcessor> processor_ptr) {
    processors_.push_back(std::move(processor_ptr));
  }

  // Parses the metrics file, passing entries to the processors that are
  // interested in them.
  void Parse();

  // Parses the metrics file, returning information about the metrics contained
  // in it.
  Manifest ParseManifest() const;

  void ClearProcessors() { processors_.clear(); }

 private:
  const std::string metrics_file_;

  // When a processor is added this class takes ownership and stores it here.
  std::vector<std::unique_ptr<MetricProcessor>> processors_;
};

// The rest of this file defines a very simple external API that can be used by
// code that only understands C. Only numeric and binary blob metrics types are
// handled. All external functions are prefixed by 'MetricsParser'.

// Helper struct, not exposed.
template <typename T>
struct ValuesAndManifest {
  using TimestampAndValue = std::pair<uint64_t, T>;

  // Values of the metric.
  std::vector<TimestampAndValue> values;
  // The manifest entry.
  PBManifestEntry manifest_entry;
};

template <typename T>
class MetricsResultHandleBase {
 public:
  // Moves the internal pointer to the next set of fields. Returns false if
  // there are no more sets of fields.
  bool Advance() {
    if (id_to_values_.empty()) {
      return false;
    } else if (next_it_ == id_to_values_.end()) {
      next_it_ = id_to_values_.begin();
      return true;
    }

    ++next_it_;
    return next_it_ != id_to_values_.end();
  }

  // Sorts all values stored by timestamp.
  void Sort() {
    for (auto& id_and_values : id_to_values_) {
      ValuesAndManifest<T>& values_and_manifest = id_and_values.second;
      auto& to_sort = values_and_manifest.values;
      std::sort(
          to_sort.begin(), to_sort.end(),
          [](const typename ValuesAndManifest<T>::TimestampAndValue& lhs,
             const typename ValuesAndManifest<T>::TimestampAndValue& rhs) {
            return lhs.first < rhs.first;
          });
    }
  }

  // Adds a new set of values to the internal map.
  void Update(uint64_t timestamp, T value, uint32_t manifest_index,
              const PBManifestEntry& manifest_entry,
              uint64_t limiting_timestamp) {
    ValuesAndManifest<T>* to_update;
    auto it = id_to_values_.find(manifest_index);
    if (it == id_to_values_.end()) {
      ValuesAndManifest<T>& new_entry = id_to_values_[manifest_index];
      new_entry.manifest_entry = manifest_entry;
      to_update = &new_entry;
    } else {
      to_update = &(it->second);
    }

    if (limiting_timestamp) {
      if (timestamp <= limiting_timestamp) {
        if (to_update->values.empty()) {
          to_update->values.emplace_back(timestamp, value);
        } else {
          typename ValuesAndManifest<T>::TimestampAndValue&
              current_timestamp_and_value = to_update->values.front();
          if (current_timestamp_and_value.first < timestamp) {
            current_timestamp_and_value.first = timestamp;
            current_timestamp_and_value.second = value;
          }
        }
      }
    } else {
      to_update->values.emplace_back(timestamp, value);
    }
  }

  // The size of (number of entries in) the current set of fields.
  uint64_t Size() {
    CHECK(next_it_ != id_to_values_.end());
    return next_it_->second.values.size();
  }

  // A string that describes the metric of the current set of fields.
  std::string MetricString() {
    CHECK(next_it_ != id_to_values_.end());
    return next_it_->second.manifest_entry.id();
  }

  // A string that describes the current set of fields.
  std::string FieldString() {
    CHECK(next_it_ != id_to_values_.end());
    return GetFieldString(next_it_->second.manifest_entry);
  }

  // Moves the values in the current set of fields to the caller.
  std::vector<typename ValuesAndManifest<T>::TimestampAndValue>&
  MutableValues() {
    CHECK(next_it_ != id_to_values_.end());
    ValuesAndManifest<T>& values_and_manifest = next_it_->second;
    return values_and_manifest.values;
  }

 protected:
  MetricsResultHandleBase() : next_it_(id_to_values_.end()) {}

  // Lists of values indexed by the manifest id.
  std::map<uint32_t, ValuesAndManifest<T>> id_to_values_;

  // Iterator to the next result. If no more results set to end.
  typename std::map<uint32_t, ValuesAndManifest<T>>::iterator next_it_;
};

// Returned from parsing functions that deal with numeric data.
class NumericMetricsResultHandle : public MetricsResultHandleBase<double> {
 public:
  // Copies all values from the current set of fields in the given arrays. They
  // should point to blocks of pre-allocated memory large enough to fit all
  // values.
  void CopyInto(uint64_t* timestamps_out, double* values_out);
};

// Returned from parsing functions that deal with generic data.
class BytesMetricsResultHandle : public MetricsResultHandleBase<std::string> {
 public:
  // Returns the size of the i-th buffer. In the numeric case (above) each value
  // is just a double, which is passed by value. In this case each value can
  // have a different size since it is an arbitrary binary blob. It is the
  // client's responsibility to ensure that the arguments for the call to
  // CopyInto are properly sized.
  uint64_t BufferSize(size_t i);

  // Copies all values from the current set of fields in the given arrays. They
  // should point to blocks of pre-allocated memory large enough to fit all
  // values. The second argument should be a buffer that can fit all bytes
  // objects from the current set of fields back to back.
  void CopyInto(uint64_t* timestamps_out, char* values_out);
};

// A convenience function that parses a metrics file and returns a map from
// <metric_id, fields_description> to a vector of <timestamp, value>. This
// function does the same as MetricsParserParse below, but is more convenient to
// use from C++ code.
std::map<std::pair<std::string, std::string>,
         std::vector<std::pair<uint64_t, double>>>
SimpleParseNumericData(const std::string& metrics_file,
                       const std::string& metric_regex,
                       const std::string& fields_to_match,
                       uint64_t min_timestamp, uint64_t max_timestamp,
                       uint64_t limiting_timestamp);

// Like above, but used when the manifest ids are known ahead of time.
std::map<std::pair<std::string, std::string>,
         std::vector<std::pair<uint64_t, double>>>
SimpleParseNumericData(const std::string& metrics_file,
                       const std::set<uint32_t> ids, uint64_t min_timestamp,
                       uint64_t max_timestamp, uint64_t limiting_timestamp);

// Will expose the parser functionality for things that are not C++
// (e.g. Python).
extern "C" {

// Returns a string summary of the manifest.
char* MetricsParserManifestSummary(const char* metrics_file);

// Returns a string summary of the values for a single metric id.
char* MetricsParserManifestMetricSummary(const char* metrics_file,
                                         const char* metric_id);

// Parses the given metrics file, looking for sets of fields that match
// 'fields_to_match' and belong to the metric(s) identified by 'metric_regex'.
// Returns a handle that can be used to read the results. The handle must be
// freed by MetricsParserResultHandleFree. Only values with timestamps more than
// or equal to min_timestamp and less than max_timestamp will be considered. If
// the limititing_timestamp argument is non-0 each set of fields will only
// contain one value -- the one that has a timestamp that is the closest to (but
// does not exceed) the limiting_timestamp argument.
NumericMetricsResultHandle* MetricsParserParse(const char* metrics_file,
                                               const char* metric_regex,
                                               const char* fields_to_match,
                                               uint64_t min_timestamp,
                                               uint64_t max_timestamp,
                                               uint64_t limiting_timestamp);

// Advances the handle to the next set of fields in the result set.
bool MetricsParserResultHandleAdvance(NumericMetricsResultHandle* handle);

// Returns the number of elements that can be extracted from the current set of
// fields.
uint64_t MetricsParserResultHandleSize(NumericMetricsResultHandle* handle);

// Returns a char* with a string describing the set of fields. The result should
// be freed with a call to MetricsParserStringFree.
char* MetricsParserResultHandleFieldString(NumericMetricsResultHandle* handle);

// Returns a char* with the id of the metric that the current set of fields
// belong to. Call MetricsParserStringFree to free.
char* MetricsParserResultHandleMetricString(NumericMetricsResultHandle* handle);

// Copies all entries associated with the current set of fields into the given
// arrays.
void MetricsParserResultHandleCopyInto(NumericMetricsResultHandle* handle,
                                       uint64_t* timestamps_out,
                                       double* values_out);

// Frees a result handle.
void MetricsParserResultHandleFree(NumericMetricsResultHandle* handle);

// Similar to MetricsParserParse, but will only handle binary blobs of data. The
// handle returned can be used to copy the binary data into externally-allocated
// buffers.
BytesMetricsResultHandle* MetricsParserBytesParse(const char* metrics_file,
                                                  const char* metric_regex,
                                                  const char* fields_to_match,
                                                  uint64_t min_timestamp,
                                                  uint64_t max_timestamp,
                                                  uint64_t limiting_timestamp);

// Advances the handle to the next set of fields in the result set.
bool MetricsParserBytesResultHandleAdvance(BytesMetricsResultHandle* handle);

// Returns the number of elements that can be extracted from the current set of
// fields.
uint64_t MetricsParserBytesResultHandleSize(BytesMetricsResultHandle* handle);

// Returns the size of the i-th field (buffer) in the current set of fields.
uint64_t MetricsParserBytesResultHandleBufferSize(
    BytesMetricsResultHandle* handle, uint64_t i);

// Copies all entries associated with the current set of fields into the given
// array back to back.
void MetricsParserBytesResultHandleCopyInto(BytesMetricsResultHandle* handle,
                                            uint64_t* timestamps_out,
                                            char* values_out);

// Returns a char* with a string describing the set of fields. The result should
// be freed with a call to MetricsParserStringFree.
char* MetricsParserBytesResultHandleFieldString(
    NumericMetricsResultHandle* handle);

// Returns a char* with the id of the metric that the current set of fields
// belong to. Call MetricsParserStringFree to free.
char* MetricsParserBytesResultHandleMetricString(
    NumericMetricsResultHandle* handle);

// Frees a result handle.
void MetricsParserBytesResultHandleFree(BytesMetricsResultHandle* handle);

// Frees a string.
void MetricsParserStringFree(char* str);

}  // extern "C"

// Like SimpleParseNumericData, but only matches data of type BinaryBlob and
// lets the user supply a function that will be used to convert the blobs.
template <typename T>
std::map<std::pair<std::string, std::string>,
         std::vector<std::pair<uint64_t, T>>>
SimpleParseCustomType(const std::string& metrics_file,
                      const std::string& metric_regex,
                      const std::string& fields_to_match,
                      uint64_t min_timestamp, uint64_t max_timestamp,
                      uint64_t limiting_timestamp,
                      std::function<T(const std::string&)> converter) {
  std::map<std::pair<std::string, std::string>,
           std::vector<std::pair<uint64_t, T>>> out;

  auto result_handle =
      std::unique_ptr<BytesMetricsResultHandle>(MetricsParserBytesParse(
          metrics_file.c_str(), metric_regex.c_str(), fields_to_match.c_str(),
          min_timestamp, max_timestamp, limiting_timestamp));
  if (!result_handle) {
    return out;
  }

  while (result_handle->Advance()) {
    std::string metric_id = result_handle->MetricString();
    std::string fields = result_handle->FieldString();
    std::vector<std::pair<uint64_t, T>>& vector = out[{metric_id, fields}];

    std::vector<std::pair<uint64_t, std::string>> vector_serialized =
        std::move(result_handle->MutableValues());
    for (const auto& timestamp_and_serialized_value : vector_serialized) {
      uint64_t timestamp = timestamp_and_serialized_value.first;
      const std::string& serialized_value =
          timestamp_and_serialized_value.second;
      vector.emplace_back(
          std::make_pair(timestamp, converter(serialized_value)));
    }
  }

  return out;
}

}  // namespace parser
}  // namsepace metrics
}  // namsepace ncode

#endif
