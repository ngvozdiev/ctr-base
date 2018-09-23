#ifndef NC_DEEP_STORAGE
#define NC_DEEP_STORAGE

#include <google/protobuf/repeated_field.h>
#include <ncode/common.h>
#include <ncode/file.h>
#include <ncode/logging.h>
#include <ncode/lru_cache.h>
#include <ncode/map_util.h>
#include <ncode/net/net_common.h>
#include <ncode/status.h>
#include <ncode/statusor.h>
#include <ncode/strutil.h>
#include <ncode/substitute.h>
#include <ncode/fwrapper.h>
#include <stddef.h>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "../../build/deep_storage.pb.h"

namespace nc {
namespace deep_storage {

class Value {
 public:
  Value(const Value& other) : type_(other.type_) {
    switch (type_) {
      case DISCRETE_NUMERIC:
        discrete_numeric_value_ = other.discrete_numeric_value_;
        break;
      case NUMERIC:
        numeric_value_ = other.numeric_value_;
        break;
      case STRING:
        new (&string_value_) std::string(other.string_value_);
        string_value_ = other.string_value_;
        break;
      case IPV4:
        ip_v4_address_ = other.ip_v4_address_;
        break;
      default:
        LOG(FATAL) << "Should not happen";
    }
  }

  ~Value() {
    if (type_ == STRING) {
      string_value_.~basic_string();
    }
  }

  Value(const char* string_value)
      : type_(STRING), string_value_(string_value) {}

  Value(const std::string& string_value)
      : type_(STRING), string_value_(string_value) {}

  Value(double numeric_value) : type_(NUMERIC), numeric_value_(numeric_value) {}

  Value(uint64_t discrete_numeric_value)
      : type_(DISCRETE_NUMERIC),
        discrete_numeric_value_(discrete_numeric_value) {}

  Value(net::IPAddress ip_address)
      : type_(IPV4), ip_v4_address_(ip_address.Raw()) {}

  Value(const ValueProto& proto) {
    switch (proto.element_case()) {
      case ValueProto::ElementCase::kDiscreteNumericValue:
        type_ = DISCRETE_NUMERIC;
        discrete_numeric_value_ = proto.discrete_numeric_value();
        break;
      case ValueProto::ElementCase::kIpV4Address:
        type_ = IPV4;
        ip_v4_address_ = net::IPAddress(proto.ip_v4_address());
        break;
      case ValueProto::ElementCase::kStringValue:
        type_ = STRING;
        string_value_ = proto.string_value();
        break;
      case ValueProto::ElementCase::kNumericValue:
        type_ = NUMERIC;
        numeric_value_ = proto.numeric_value();
        break;
      case ValueProto::ELEMENT_NOT_SET:
        LOG(FATAL) << "No element set";
    }
  }

  std::string ToString() const {
    switch (type_) {
      case DISCRETE_NUMERIC:
        return std::to_string(discrete_numeric_value_);
      case NUMERIC:
        return std::to_string(numeric_value_);
      case STRING:
        return string_value_;
      case IPV4:
        return net::IPToStringOrDie(ip_v4_address_);
      default:
        LOG(FATAL) << "Should not happen";
        return "";
    }
  }

  ValueProto ToProto() const {
    ValueProto out;

    switch (type_) {
      case DISCRETE_NUMERIC:
        out.set_discrete_numeric_value(discrete_numeric_value_);
        break;
      case NUMERIC:
        out.set_numeric_value(numeric_value_);
        break;
      case STRING:
        out.set_string_value(string_value_);
        break;
      case IPV4:
        out.set_ip_v4_address(ip_v4_address_.Raw());
        break;
      default:
        LOG(FATAL) << "Should not happen";
    }

    return out;
  }

 private:
  enum { DISCRETE_NUMERIC, NUMERIC, STRING, IPV4 } type_;

  union {
    uint64_t discrete_numeric_value_;
    double numeric_value_;
    std::string string_value_;
    net::IPAddress ip_v4_address_;
  };
};

// An entry is just a byte vector.
using Entry = std::vector<uint8_t>;

// Set of keys.
using KeySet = std::set<std::string>;

// A map from keys to values.
using KeyValueMap = std::map<std::string, Value>;

class EntryHandle {
 public:
  EntryHandle(const std::vector<Value>& values,
              std::chrono::microseconds timestamp, uint64_t chunk_id,
              uint64_t offset_in_chunk, uint64_t entry_size)
      : values_(values),
        timestamp_(timestamp),
        chunk_id_(chunk_id),
        offset_in_chunk_(offset_in_chunk),
        entry_size_(entry_size) {}

  EntryHandle(const EntryHandleProto& proto, uint64_t chunk_id)
      : timestamp_(std::chrono::microseconds(proto.creation_timestamp())),
        chunk_id_(chunk_id),
        offset_in_chunk_(proto.offset_in_chunk()),
        entry_size_(proto.entry_size()) {
    for (const auto& v : proto.values()) {
      values_.emplace_back(v);
    }
  }

  EntryHandleProto ToProto() const {
    EntryHandleProto out;
    out.set_creation_timestamp(timestamp_.count());
    out.set_offset_in_chunk(offset_in_chunk_);
    out.set_entry_size(entry_size_);
    for (const auto& v : values_) {
      *out.add_values() = v.ToProto();
    }

    return out;
  }

  uint64_t offset_in_chunk() const { return offset_in_chunk_; }

  uint64_t entry_size() const { return entry_size_; }

  const std::vector<Value>& values() { return values_; }

  uint64_t chunk_id() const { return chunk_id_; }

 private:
  // The values for this entry.
  std::vector<Value> values_;

  // Time the entry was created.
  std::chrono::microseconds timestamp_;

  // Identifies the chunk.
  uint64_t chunk_id_;

  // Identifies the entry in the chunk.
  uint64_t offset_in_chunk_;
  uint64_t entry_size_;

  DISALLOW_COPY_AND_ASSIGN(EntryHandle);
};

class FileOpener {
 public:
  FileOpener(uint32_t max_open_files, const std::string& pattern)
      : pattern_(pattern), open_files_(max_open_files) {
    CHECK(max_open_files > 0);
  }

  nc::StatusOr<FWrapper*> GetStream(uint64_t id);

 private:
  const std::string pattern_;
  nc::LRUCache<uint64_t, FWrapper> open_files_;
};

class Chunk {
 public:
  static nc::StatusOr<std::unique_ptr<Chunk>> ReadChunk(FileOpener* file_opener,
                                                        uint64_t chunk_id);

  static nc::StatusOr<std::unique_ptr<Chunk>> InitChunk(FileOpener* file_opener,
                                                        const KeySet& keys,
                                                        uint64_t chunk_id);

  Chunk(FileOpener* file_opener, uint64_t chunk_id, const KeySet& keys,
        std::vector<std::unique_ptr<EntryHandle>>* entries = nullptr)
      : chunk_id_(chunk_id), keys_(keys), file_opener_(file_opener) {
    if (entries != nullptr) {
      entries_ = std::move(*entries);
    }
  }

  ~Chunk() {
    auto result = file_opener_->GetStream(chunk_id_);
    if (result.ok()) {
      result.ValueOrDie()->Flush();
    }
  }

  // Appends an entry and returns its handle.
  StatusOr<const EntryHandle*> AppendEntry(const std::vector<Value>& values,
                                           const Entry& entry);

  StatusOr<Entry> ReadEntry(const EntryHandle& handle);

  const KeySet& keys() const { return keys_; }

 private:
  // Unique identifier for the chunk.
  const uint64_t chunk_id_;

  // Keys for all values in the chunk.
  const KeySet keys_;

  // Entries in the chunk.
  std::vector<std::unique_ptr<EntryHandle>> entries_;

  // Read/append files come from here.
  FileOpener* file_opener_;

  DISALLOW_COPY_AND_ASSIGN(Chunk);
};

class Storage {
 public:
  Storage(const std::string& chunk_pattern, uint32_t max_open_files = 100)
      : file_opener_(max_open_files, chunk_pattern) {}

  StatusOr<const EntryHandle*> Put(const KeyValueMap& kv_map,
                                   const Entry& entry) {
    KeySet keys;
    std::vector<Value> values;
    for (const auto& key_and_value : kv_map) {
      keys.insert(key_and_value.first);
      values.emplace_back(key_and_value.second);
    }

    Chunk* chunk;
    ASSIGN_OR_RETURN(chunk, GetChunk(keys));
    return chunk->AppendEntry(values, entry);
  }

  StatusOr<Entry> Get(const EntryHandle& handle) {
    Chunk* chunk = FindSmartPtrOrNull(chunks_, handle.chunk_id());
    if (chunk == nullptr) {
      return Status(error::INVALID_ARGUMENT,
                    StrCat("Bad chunk id ", handle.chunk_id()));
    }

    return chunk->ReadEntry(handle);
  }

 private:
  static constexpr char kFilePattern[] = "chunk_$0";
  static constexpr uint32_t kMaxOpenFiles = 100;

  uint64_t PickFreeChunkId() {
    return GetRandomId([this](uint64_t candidate) {
      return File::Exists(Substitute(kFilePattern, candidate));
    }, &rnd_);
  }

  StatusOr<Chunk*> GetChunk(const KeySet& keys) {
    Chunk* chunks_for_keys = FindPtrOrNull(keys_to_chunks_, keys);
    if (chunks_for_keys != nullptr) {
      return chunks_for_keys;
    }

    // Unable to find existing chunks, will allocate a new one.
    uint64_t new_id = PickFreeChunkId();

    auto init_result = Chunk::InitChunk(&file_opener_, keys, new_id);
    if (!init_result.ok()) {
      return init_result.status();
    }

    std::unique_ptr<Chunk> new_chunk = init_result.ConsumeValueOrDie();
    Chunk* chunk_ptr = new_chunk.get();
    chunks_[new_id] = std::move(new_chunk);
    keys_to_chunks_[keys] = chunk_ptr;
    return chunk_ptr;
  }

  // Read/append files come from here.
  FileOpener file_opener_;

  // The chunks.
  std::map<uint64_t, std::unique_ptr<Chunk>> chunks_;

  // A map from a set of keys to all chunks that hold values for those keys.
  std::map<KeySet, Chunk*> keys_to_chunks_;
  std::mt19937 rnd_;
};

}  // namespace deep_storage
}  // namespace nc

#endif
