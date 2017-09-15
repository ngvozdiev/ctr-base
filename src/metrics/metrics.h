#ifndef NCODE_METRICS_METRIC_H
#define NCODE_METRICS_METRIC_H

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/repeated_field.h>
#include <stddef.h>
#include <stdint.h>
#include <cassert>
#include <chrono>
#include <cstdbool>
#include <ctime>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <mutex>
#include <atomic>

#include "metrics.pb.h"
#include "ncode_common/src/circular_array.h"
#include "ncode_common/src/common.h"
#include "ncode_common/src/event_queue.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/stats.h"
#include "ncode_common/src/strutil.h"

namespace nc {
namespace metrics {

// A file descriptor along with an output stream. All are populated upon
// construction.
class OutputStream {
 public:
  OutputStream(const std::string& file);

  ~OutputStream();

  // Writes a series of entries to the stream.
  void WriteBulk(const std::vector<PBMetricEntry>& entries,
                 uint32_t manifest_index);

  // Writes a single entry to the stream.
  void WriteSingle(const PBMetricEntry& entry, uint32_t manifest_index);

 private:
  // Writes a protobuf to the stream.
  bool WriteDelimitedTo(const PBMetricEntry& entry, uint32_t manifest_index);

  // A file descriptor. Closed on destruction.
  int fd_;

  // File output stream. Owned by this object.
  std::unique_ptr<google::protobuf::io::FileOutputStream> file_output_;

  // The output stream is thread safe.
  std::mutex mu_;

  DISALLOW_COPY_AND_ASSIGN(OutputStream);
};

template <typename T>
struct Entry {
  T value;             // The actual value.
  uint64_t timestamp;  // Creation time.

  bool operator==(const Entry& rhs) const {
    return (value == rhs.value) && (timestamp == rhs.timestamp);
  }

  bool operator!=(const Entry& rhs) const { return !operator==(rhs); }
};

template <typename T>
void SaveEntryToProtobuf(const Entry<T>& entry, PBMetricEntry* out) {
  Unused(entry);
  Unused(out);
  assert(false);
}

template <>
void SaveEntryToProtobuf<uint64_t>(const Entry<uint64_t>& entry,
                                   PBMetricEntry* out);

template <>
void SaveEntryToProtobuf<uint32_t>(const Entry<uint32_t>& entry,
                                   PBMetricEntry* out);

template <>
void SaveEntryToProtobuf<bool>(const Entry<bool>& entry, PBMetricEntry* out);

template <>
void SaveEntryToProtobuf<double>(const Entry<double>& entry,
                                 PBMetricEntry* out);

template <>
void SaveEntryToProtobuf<std::string>(const Entry<std::string>& entry,
                                      PBMetricEntry* out);

template <>
void SaveEntryToProtobuf<BytesBlob>(const Entry<BytesBlob>& entry,
                                    PBMetricEntry* out);

template <>
void SaveEntryToProtobuf<nc::DiscreteDistribution<uint64_t>>(
    const Entry<nc::DiscreteDistribution<uint64_t>>& entry, PBMetricEntry* out);

template <>
void SaveEntryToProtobuf<nc::DiscreteDistribution<int64_t>>(
    const Entry<nc::DiscreteDistribution<int64_t>>& entry, PBMetricEntry* out);

// A generic interface for a class that knows how to provide timestamps.
class TimestampProviderInterface {
 public:
  virtual ~TimestampProviderInterface() {}
  virtual uint64_t GetTimestamp() const = 0;
  virtual const char* TimestampUnits() const = 0;
  virtual std::string TimestampToString(uint64_t timestamp) const = 0;
};

class DefaultTimestampProvider : public TimestampProviderInterface {
 public:
  static constexpr const char* kNanoseconds = "nanoseconds";
  DefaultTimestampProvider() {}

  uint64_t GetTimestamp() const override {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    auto now = high_resolution_clock::now();
    return duration_cast<std::chrono::nanoseconds>(now.time_since_epoch())
        .count();
  }

  const char* TimestampUnits() const override { return kNanoseconds; }

  std::string TimestampToString(uint64_t timestamp) const override {
    using std::chrono::system_clock;
    const auto dt = std::chrono::nanoseconds(timestamp);
    const std::chrono::time_point<system_clock> tp_after_duration(
        std::chrono::duration_cast<std::chrono::milliseconds>(dt));
    time_t time_after_duration = system_clock::to_time_t(tp_after_duration);
    uint64_t milliseconds_remainder = (timestamp / 1000) % 1000;

    char s[128];
    strftime(s, 128, "%H:%M:%S ", std::localtime(&time_after_duration));

    std::string time_str(s);
    time_str += std::to_string(milliseconds_remainder) + "ms";
    return time_str;
  }
};

class MetricBase;

class MetricHandleBase {
 public:
  // Two timestamps.
  using Span = std::pair<uint64_t, uint64_t>;

  virtual ~MetricHandleBase() {}

  // The application is polled for the most recent value and it is added to the
  // history. If there is no application callback this does nothing.
  virtual void Poll() = 0;

  // Causes in-mem values to get flushed to the persistence file. If no
  // persistence is enabled does not have effect.
  virtual void Persist() = 0;

  // Whether or not the metric has a callback that will produce values upon
  // polling.
  bool has_callback() const { return has_callback_; }

  // Number of elements in the history.
  virtual size_t NumEntriesInMem() const = 0;

  // The timestamps of the first and the last element of this metric's history.
  virtual Span EntriesInMemSpan() const = 0;

  size_t metric_index() const { return metric_index_; }

  std::string TimestampToString(uint64_t timestamp) const;

  void set_metric_index(size_t index) { metric_index_ = index; }

 protected:
  MetricHandleBase(bool has_callback, MetricBase* parent_metric);

  // Is there a callback in the implementing class.
  const bool has_callback_;

  // This metric's index.
  size_t metric_index_;

  MetricBase* parent_metric_;

 private:
  DISALLOW_COPY_AND_ASSIGN(MetricHandleBase);
};

// Mixin for classes that need to have a mutex.
template <bool>
struct SomethingWithMutex {
  inline std::unique_lock<std::mutex> GetLock() const {
    return std::unique_lock<std::mutex>();
  }
};
template <>
struct SomethingWithMutex<true> {
  inline std::unique_lock<std::mutex> GetLock() const {
    return std::unique_lock<std::mutex>(mu_);
  }

  mutable std::mutex mu_;
};

template <typename EntryType, bool ThreadSafe>
class MetricHandle : public MetricHandleBase,
                     private SomethingWithMutex<ThreadSafe> {
 public:
  static constexpr size_t kEntriesInMem = 32;
  using MetricCallback = std::function<EntryType()>;

  MetricHandle(MetricBase* parent_metric, MetricCallback callback);

  MetricHandle(MetricBase* parent_metric);

  size_t NumEntriesInMem() const override {
    std::unique_lock<std::mutex> lock = GetLock();
    return storage_.size();
  }

  Span EntriesInMemSpan() const override;

  // Explicitly adds a value to the metric. The timestamp associated with this
  // value will be taken from the MetricProvider's timestamp generator.
  Entry<EntryType> AddValue(EntryType value);

  // Like AddValue, but the timestamp can be specified.
  Entry<EntryType> AddValueWithTimestamp(uint64_t time_now, EntryType value);

  void Poll() override;

  void Persist() override;

  bool MostRecentEntryInMem(Entry<EntryType>* entry) const {
    std::unique_lock<std::mutex> lock = GetLock();
    if (storage_.empty()) {
      return false;
    }

    *entry = storage_.MostRecentValueOrDie();
    return true;
  }

 private:
  using SomethingWithMutex<ThreadSafe>::GetLock;

  std::vector<Entry<EntryType>> PersistPrivate();

  void ValuesToDisk(const std::vector<Entry<EntryType>>& values);

  Entry<EntryType> PollAndAdd(uint64_t time_now);

  Entry<EntryType> AddValuePrivate(EntryType value, uint64_t time_now);

  // A callback to the app to provide values.
  const MetricCallback app_callback_;

  // If logging is enabled this counter represents the number of log records
  // that can be skipped before the log is flushed. Every time a new message is
  // logged this is decremented and when it reaches 0 it the log is flushed. The
  // counter is reset to the max number of messages that can be stored in the
  // in-mem circular buffer (NumValues).
  size_t log_tokens_;

  // Elements stored.
  CircularArray<Entry<EntryType>, kEntriesInMem> storage_;

  DISALLOW_COPY_AND_ASSIGN(MetricHandle);
};

template <typename T>
using UnsafeMetricHandle = MetricHandle<T, false>;

void PopulateManifestEntryType(PBManifestEntry* out, uint64_t* dummy);
void PopulateManifestEntryType(PBManifestEntry* out, uint32_t* dummy);
void PopulateManifestEntryType(PBManifestEntry* out, bool* dummy);
void PopulateManifestEntryType(PBManifestEntry* out, std::string* dummy);
void PopulateManifestEntryType(PBManifestEntry* out, double* dummy);
void PopulateManifestEntryType(PBManifestEntry* out, BytesBlob* dummy);
void PopulateManifestEntryType(PBManifestEntry* out,
                               nc::DiscreteDistribution<uint64_t>* dummy);
void PopulateManifestEntryType(PBManifestEntry* out,
                               nc::DiscreteDistribution<int64_t>* dummy);

void PopulateManifestEntryField(PBMetricField* field, uint64_t value);
void PopulateManifestEntryField(PBMetricField* field, uint32_t value);
void PopulateManifestEntryField(PBMetricField* field, bool value);
void PopulateManifestEntryField(PBMetricField* field, const std::string& value);

template <size_t remaining, typename... FieldTypes>
struct PopulateManifestEntry {
  static void Execute(PBManifestEntry* out,
                      const std::tuple<FieldTypes...>& fields) {
    constexpr size_t num_fields =
        std::tuple_size<std::tuple<FieldTypes...>>::value;
    PBMetricField* new_field = out->mutable_fields(num_fields - remaining);
    const auto& field_value = std::get<num_fields - remaining>(fields);
    PopulateManifestEntryField(new_field, field_value);
    PopulateManifestEntry<remaining - 1, FieldTypes...>::Execute(out, fields);
  }
};

template <typename... FieldTypes>
struct PopulateManifestEntry<0, FieldTypes...> {
  static void Execute(PBManifestEntry* out,
                      const std::tuple<FieldTypes...>& fields) {
    Unused(out);
    Unused(fields);
  }
};

class MetricManager;

// All metrics inherit from this class.
class MetricBase {
 public:
  static constexpr uint32_t kManifestEntryMetaIndex =
      std::numeric_limits<uint32_t>::max();

  virtual ~MetricBase() {}

  // Returns PBManifest entries for all handles grouped by index in the
  // manifest.
  virtual std::map<size_t, std::unique_ptr<PBManifestEntry>>
  ManifestIndexToManifestEntry() const = 0;

  // Causes all handles to be persisted.
  virtual void PersistAllHandles() = 0;

  // Causes all handles that have callbacks registered to get their value.
  virtual void PollAllHandles() = 0;

  // Clears all handles.
  virtual void ClearAllHandles() = 0;

  void SetLocalOutputStream(std::unique_ptr<OutputStream> output_stream);

  std::string id() const { return base_entry_.id(); }

  // Returns the local output stream if not null, or the parent's input stream
  // or null if no output stream is set. After calling this function the stream
  // cannot be set.
  OutputStream* OutputStreamOrNull();

  const TimestampProviderInterface* timestamp_provider() const;

  // Whether or not the output stream can be set. Attempting to set the output
  // stream when it cannot be set will result in a crash.
  bool stream_locked() const { return stream_locked_; }

 protected:
  MetricBase(MetricManager* metric_manager, PBManifestEntry base_entry)
      : parent_manager_(metric_manager),
        base_entry_(base_entry),
        local_current_index_(std::numeric_limits<size_t>::max()),
        stream_locked_(false) {}

  size_t NextIndex();

  // The parent manager.
  MetricManager* const parent_manager_;

  // All MetricHandles that are produced by the metric will have a copy of this
  // PBManifestEntry, but with different values.
  PBManifestEntry base_entry_;

 private:
  // Only set if the metric has its own output file. If it does not it will use
  // the parent manager's output stream.
  std::unique_ptr<OutputStream> local_output_stream_;

  // Only used if 'local_output_stream' is not null.
  size_t local_current_index_;

  // Set to true if no more changes to the output stream are allowed.
  std::atomic<bool> stream_locked_;
};

// A metric is type associated with a combination of fields.
template <typename EntryType, bool ThreadSafe, typename... FieldTypes>
class Metric : public MetricBase, private SomethingWithMutex<ThreadSafe> {
 public:
  using MetricCallback = std::function<EntryType()>;
  using HandleType = MetricHandle<EntryType, ThreadSafe>;

  Metric(MetricManager* metric_manager, PBManifestEntry base_entry)
      : MetricBase(metric_manager, base_entry) {}

  void PersistAllHandles() override {
    std::unique_lock<std::mutex> lock = GetLock();
    for (auto& handle_fields_and_handle : fields_to_handle_) {
      HandleType& handle = handle_fields_and_handle.second;
      handle.Persist();
    }
  }

  void PollAllHandles() override {
    std::unique_lock<std::mutex> lock = GetLock();
    for (auto& handle_fields_and_handle : fields_to_handle_) {
      HandleType& handle = handle_fields_and_handle.second;
      handle.Poll();
    }
  }

  void ClearAllHandles() override {
    std::unique_lock<std::mutex> lock = GetLock();
    fields_to_handle_.clear();
  }

  // Gets a handle that can be used to add values to this metric. The handle
  // pointer is owned by this class.
  HandleType* GetHandle(FieldTypes... fields);

  HandleType* GetHandle(MetricCallback callback, FieldTypes... fields);

  // Returns PBManifest entries for all handles grouped by index in the
  // manifest.
  std::map<size_t, std::unique_ptr<PBManifestEntry>>
  ManifestIndexToManifestEntry() const override {
    std::unique_lock<std::mutex> lock = GetLock();
    std::map<size_t, std::unique_ptr<PBManifestEntry>> return_map;
    for (const auto& fields_and_handle : fields_to_handle_) {
      const std::tuple<FieldTypes...>& fields = fields_and_handle.first;
      const HandleType& handle = fields_and_handle.second;

      // Each handle will get a copy of the base manifest entry, but modified
      // for its fields.
      auto entry_for_handle = make_unique<PBManifestEntry>(base_entry_);
      PopulateManifestEntry<kNumFields, FieldTypes...>::Execute(
          entry_for_handle.get(), fields);
      return_map[handle.metric_index()] = std::move(entry_for_handle);
    }

    return return_map;
  }

 private:
  using SomethingWithMutex<ThreadSafe>::GetLock;

  static constexpr size_t kNumFields =
      std::tuple_size<std::tuple<FieldTypes...>>::value;

  std::map<std::tuple<FieldTypes...>, HandleType> fields_to_handle_;

  DISALLOW_COPY_AND_ASSIGN(Metric);
};

class MetricManager {
 public:
  MetricManager();

  ~MetricManager();

  void PersistAllMetrics();

  void PollAllMetrics();

  // Sets the output. Should only be called once. Each metric will have a
  // separate file, named after the metric id. Those files will be stored in a
  // directory 'output'.
  void SetOutput(const std::string& output);

  template <typename EntryType, bool ThreadSafe, typename... FieldTypes,
            typename... DescTypes>
  Metric<EntryType, ThreadSafe, FieldTypes...>* GetMetric(
      const std::string& id, const std::string& metric_description,
      DescTypes... field_descriptions) {
    static_assert(sizeof...(FieldTypes) == sizeof...(DescTypes),
                  "Fields / field descriptions count mismatch");
    PBManifestEntry new_manifest_base_entry;
    new_manifest_base_entry.set_id(id);
    new_manifest_base_entry.set_description(metric_description);
    PopulateManifestEntryType(&new_manifest_base_entry,
                              static_cast<EntryType*>(0));
    PopulateManifestEntryPrivate(&new_manifest_base_entry,
                                 field_descriptions...);

    auto metric = make_unique<Metric<EntryType, ThreadSafe, FieldTypes...>>(
        this, new_manifest_base_entry);
    auto* raw_metric_ptr = metric.get();

    std::lock_guard<std::mutex> lock(mu_);
    if (!output_directory_.empty()) {
      std::string local_output = StrCat(output_directory_, "/", id);
      auto local_output_stream = make_unique<OutputStream>(local_output);
      metric->SetLocalOutputStream(std::move(local_output_stream));
    }

    all_metrics_.push_back(std::move(metric));
    return raw_metric_ptr;
  }

  template <typename EntryType, typename... FieldTypes, typename... DescTypes>
  Metric<EntryType, true, FieldTypes...>* GetThreadSafeMetric(
      const std::string& id, const std::string& metric_description,
      DescTypes... field_descriptions) {
    return GetMetric<EntryType, true, FieldTypes...>(id, metric_description,
                                                     field_descriptions...);
  }

  template <typename EntryType, typename... FieldTypes, typename... DescTypes>
  Metric<EntryType, false, FieldTypes...>* GetUnsafeMetric(
      const std::string& id, const std::string& metric_description,
      DescTypes... field_descriptions) {
    return GetMetric<EntryType, false, FieldTypes...>(id, metric_description,
                                                      field_descriptions...);
  }

  const TimestampProviderInterface* timestamp_provider() const;

  void Clear() {
    std::lock_guard<std::mutex> lock(mu_);
    PersistAllMetrics();
    timestamp_provider_ = make_unique<DefaultTimestampProvider>();
    for (const auto& metric : all_metrics_) {
      metric->ClearAllHandles();
    }
  }

  void set_timestamp_provider(
      std::unique_ptr<TimestampProviderInterface> timestamp_provider) {
    Clear();
    timestamp_provider_ = std::move(timestamp_provider);
  }

 private:
  template <typename FieldType, typename... FieldTypes>
  void PopulateManifestEntryPrivate(PBManifestEntry* out,
                                    FieldType next_field_description,
                                    FieldTypes... field_descriptions) {
    PBMetricField* new_field = out->add_fields();
    new_field->set_description(next_field_description);

    return PopulateManifestEntryPrivate(out, field_descriptions...);
  }

  void PopulateManifestEntryPrivate(PBManifestEntry* out) {
    Unused(out);
    return;
  }

  // A class that knows how to provide timestamps.
  std::unique_ptr<TimestampProviderInterface> timestamp_provider_;

  // All metrics are stored here.
  std::vector<std::unique_ptr<MetricBase>> all_metrics_;

  // Output directory where the per-metric logs are kept. Only non-empty if
  // there are per-metric logs (output_stream_ is null).
  std::string output_directory_;

  // Protects current_index_ and all_metrics_;
  std::mutex mu_;
};

// Returns the default metric manager singleton.
MetricManager* DefaultMetricManager();

// Should be called before any calls to DefaultMetricManager.
void InitMetrics();

template <typename EntryType, bool ThreadSafe>
MetricHandle<EntryType, ThreadSafe>::MetricHandle(
    MetricBase* parent_metric, std::function<EntryType()> callback)
    : MetricHandleBase(true, parent_metric),
      app_callback_(callback),
      log_tokens_(kEntriesInMem) {
  CHECK(callback.operator bool() == true);
}

template <typename EntryType, bool ThreadSafe>
MetricHandle<EntryType, ThreadSafe>::MetricHandle(MetricBase* parent_metric)
    : MetricHandleBase(false, parent_metric), log_tokens_(kEntriesInMem) {}

template <typename EntryType, bool ThreadSafe>
Entry<EntryType> MetricHandle<EntryType, ThreadSafe>::AddValue(
    EntryType value) {
  uint64_t time_now = parent_metric_->timestamp_provider()->GetTimestamp();
  return AddValuePrivate(value, time_now);
}

template <typename EntryType, bool ThreadSafe>
Entry<EntryType> MetricHandle<EntryType, ThreadSafe>::AddValueWithTimestamp(
    uint64_t time_now, EntryType value) {
  return AddValuePrivate(value, time_now);
}

template <typename EntryType, bool ThreadSafe>
void MetricHandle<EntryType, ThreadSafe>::Poll() {
  if (!app_callback_) {
    return;
  }

  uint64_t time_now = parent_metric_->timestamp_provider()->GetTimestamp();
  PollAndAdd(time_now);
}

template <typename EntryType, bool ThreadSafe>
typename MetricHandleBase::Span
MetricHandle<EntryType, ThreadSafe>::EntriesInMemSpan() const {
  std::unique_lock<std::mutex> lock = GetLock();
  if (storage_.empty()) {
    return std::make_pair(std::numeric_limits<uint64_t>::min(),
                          std::numeric_limits<uint64_t>::min());
  }

  const Entry<EntryType>& min_value = storage_.OldestValueOrDie();
  const Entry<EntryType>& max_value = storage_.MostRecentValueOrDie();
  return std::make_pair(min_value.timestamp, max_value.timestamp);
}

template <typename EntryType, bool ThreadSafe>
void MetricHandle<EntryType, ThreadSafe>::Persist() {
  std::vector<Entry<EntryType>> values = PersistPrivate();
  ValuesToDisk(values);
}

template <typename EntryType, bool ThreadSafe>
std::vector<Entry<EntryType>>
MetricHandle<EntryType, ThreadSafe>::PersistPrivate() {
  std::vector<Entry<EntryType>> values = storage_.GetValues();

  // How many values we have to skip from the start of the log.
  size_t values_to_log = kEntriesInMem - log_tokens_;
  DCHECK(values_to_log == values.size());
  log_tokens_ = kEntriesInMem;
  return values;
}

template <typename EntryType, bool ThreadSafe>
void MetricHandle<EntryType, ThreadSafe>::ValuesToDisk(
    const std::vector<Entry<EntryType>>& values) {
  OutputStream* output_stream = parent_metric_->OutputStreamOrNull();
  if (!output_stream) {
    return;
  }

  std::vector<PBMetricEntry> entries_to_stream(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const Entry<EntryType>& value = values[i];
    SaveEntryToProtobuf(value, &entries_to_stream[i]);
  }
  output_stream->WriteBulk(entries_to_stream, metric_index_);
}

template <typename EntryType, bool ThreadSafe>
Entry<EntryType> MetricHandle<EntryType, ThreadSafe>::PollAndAdd(
    uint64_t time_now) {
  EntryType value = app_callback_();
  Entry<EntryType> return_value = AddValuePrivate(value, time_now);
  return return_value;
}

template <typename EntryType, bool ThreadSafe>
Entry<EntryType> MetricHandle<EntryType, ThreadSafe>::AddValuePrivate(
    EntryType value, uint64_t time_now) {
  Entry<EntryType> entry = {value, time_now};
  std::unique_lock<std::mutex> lock = GetLock();
  storage_.AddValue(entry);

  --log_tokens_;
  if (log_tokens_ == 0) {
    Persist();
  }

  return entry;
}

template <>
Entry<double> MetricHandle<double, true>::AddValuePrivate(double value,
                                                          uint64_t time_now);

template <>
Entry<double> MetricHandle<double, false>::AddValuePrivate(double value,
                                                           uint64_t time_now);

template <typename EntryType, bool ThreadSafe, typename... FieldTypes>
typename Metric<EntryType, ThreadSafe, FieldTypes...>::HandleType* Metric<
    EntryType, ThreadSafe, FieldTypes...>::GetHandle(FieldTypes... fields) {
  std::unique_lock<std::mutex> lock = GetLock();
  OutputStream* output_stream = OutputStreamOrNull();

  auto it = fields_to_handle_.emplace(std::piecewise_construct,
                                      std::forward_as_tuple(fields...),
                                      std::forward_as_tuple(this));
  HandleType* handle = &(it.first->second);
  if (it.second) {
    handle->set_metric_index(NextIndex());

    PBMetricEntry entry;
    *entry.mutable_manifest_entry() = base_entry_;
    PopulateManifestEntry<kNumFields, FieldTypes...>::Execute(
        entry.mutable_manifest_entry(), std::forward_as_tuple(fields...));
    if (output_stream) {
      output_stream->WriteSingle(entry, kManifestEntryMetaIndex);
    }
  }

  return handle;
}

template <typename EntryType, bool ThreadSafe, typename... FieldTypes>
typename Metric<EntryType, ThreadSafe, FieldTypes...>::HandleType*
Metric<EntryType, ThreadSafe, FieldTypes...>::GetHandle(MetricCallback callback,
                                                        FieldTypes... fields) {
  std::unique_lock<std::mutex> lock = GetLock();
  OutputStream* output_stream = OutputStreamOrNull();

  auto it = fields_to_handle_.emplace(std::piecewise_construct,
                                      std::forward_as_tuple(fields...),
                                      std::forward_as_tuple(this, callback));
  HandleType* handle = &(it.first->second);
  if (it.second) {
    size_t next_index = NextIndex();
    handle->set_metric_index(next_index);

    PBMetricEntry entry;
    *entry.mutable_manifest_entry() = base_entry_;
    PopulateManifestEntry<kNumFields, FieldTypes...>::Execute(
        entry.mutable_manifest_entry(), std::forward_as_tuple(fields...));
    if (output_stream) {
      output_stream->WriteSingle(entry, kManifestEntryMetaIndex);
    }
  }
  return handle;
}

// Polls all metrics in the default metric manager regularly.
class DefaultMetricManagerPoller : public EventConsumer {
 public:
  DefaultMetricManagerPoller(std::chrono::milliseconds period,
                             EventQueue* event_queue);

  void HandleEvent() override;

  void EnqueueNext() { EnqueueIn(event_queue()->ToTime(period_)); }

 private:
  std::chrono::milliseconds period_;
};

// A timestamp provider that always returns 0, effectively switching off
// timestamps.
class NullTimestampProvider : public TimestampProviderInterface {
  static constexpr const char* kPicoseconds = "picoseconds";

  uint64_t GetTimestamp() const override { return 0; }

  const char* TimestampUnits() const override { return kPicoseconds; }

  std::string TimestampToString(uint64_t timestamp) const override {
    CHECK(timestamp == 0);
    return "0";
  }
};

// A timestamp provider that is based off an event queue.
class SimTimestampProvider : public TimestampProviderInterface {
 public:
  static constexpr const char* kPicoseconds = "picoseconds";
  SimTimestampProvider(SimTimeEventQueue* event_queue)
      : event_queue_(event_queue) {}

  uint64_t GetTimestamp() const override {
    return event_queue_->CurrentTime().Raw();
  }

  const char* TimestampUnits() const override { return kPicoseconds; }

  std::string TimestampToString(uint64_t timestamp) const override {
    using std::chrono::system_clock;
    const auto dt = SimTimeEventQueue::picoseconds(timestamp);
    const std::chrono::time_point<system_clock> tp_after_duration(
        std::chrono::duration_cast<std::chrono::milliseconds>(dt));
    time_t time_after_duration = system_clock::to_time_t(tp_after_duration);
    uint64_t milliseconds_remainder = (timestamp / 1000) % 1000;

    char s[128];
    strftime(s, 128, "%H:%M:%S ", std::localtime(&time_after_duration));

    std::string time_str(s);
    time_str += std::to_string(milliseconds_remainder) + "ms";
    return time_str;
  }

 private:
  SimTimeEventQueue* event_queue_;
};

}  // namespace metrics
}  // namespace nc

#endif /* NCODE_METRICS_METRIC_H */
