#include "deep_storage.h"

#include <ncode/file.h>
#include <ncode/port.h>
#include <ncode/status.h>
#include <ncode/strutil.h>
#include <sys/errno.h>
#include <chrono>
#include <cstring>
#include <string>

namespace nc {
namespace deep_storage {

static std::chrono::microseconds TimeNow() {
  return std::chrono::duration_cast<std::chrono::microseconds>(
      std::chrono::high_resolution_clock::now().time_since_epoch());
}

StatusOr<FWrapper*> FileOpener::GetStream(uint64_t id) {
  FWrapper* out = open_files_.FindOrNull(id);
  if (out != nullptr) {
    return out;
  }

  std::string name = Substitute(pattern_.c_str(), id);
  StatusOr<FWrapper> status_or_wrapper = FWrapper::Open(name, "a+");
  if (!status_or_wrapper.ok()) {
    return status_or_wrapper.status();
  }

  std::unique_ptr<FWrapper> new_wrapper =
      make_unique<FWrapper>(status_or_wrapper.ConsumeValueOrDie());
  FWrapper& fw = open_files_.InsertNew(id, std::move(new_wrapper));
  return &fw;
}

template <typename T>
static Status ReadProto(FWrapper* file_w, size_t message_size, T* out) {
  std::string data;
  data.resize(message_size);

  RETURN_IF_ERROR(file_w->Read(message_size, &(data[0])));
  out->Clear();
  if (!out->ParseFromString(data)) {
    return Status(error::INTERNAL, StrCat("Unable to deserialize"));
  }

  return Status::OK;
}

static StatusOr<std::unique_ptr<EntryHandle>> ReadChunkStatic(
    uint64_t chunk_id, FWrapper* file_w, EntryHandleProto* tmp) {
  uint32_t entry_handle_size;
  ASSIGN_OR_RETURN(entry_handle_size, file_w->ReadUint32());

  RETURN_IF_ERROR(ReadProto<EntryHandleProto>(file_w, entry_handle_size, tmp));
  RETURN_IF_ERROR(file_w->Skip(tmp->entry_size()));
  return make_unique<EntryHandle>(*tmp, chunk_id);
}

StatusOr<std::unique_ptr<Chunk>> Chunk::ReadChunk(FileOpener* file_opener,
                                                  uint64_t chunk_id) {
  Timer timer;

  FWrapper* file_w;
  ASSIGN_OR_RETURN(file_w, file_opener->GetStream(chunk_id));
  RETURN_IF_ERROR(file_w->Seek(0));

  uint32_t header_size;
  ASSIGN_OR_RETURN(header_size, file_w->ReadUint32());

  ChunkFileHeaderProto header;
  RETURN_IF_ERROR(
      ReadProto<ChunkFileHeaderProto>(file_w, header_size, &header));

  std::vector<std::unique_ptr<EntryHandle>> entries;
  EntryHandleProto tmp;
  while (true) {
    auto read_output = ReadChunkStatic(chunk_id, file_w, &tmp);
    if (!read_output.ok()) {
      break;
    }

    entries.emplace_back(read_output.ConsumeValueOrDie());
  }

  uint64_t bytes_read;
  ASSIGN_OR_RETURN(bytes_read, file_w->Tell());

  uint64_t total_file_size;
  ASSIGN_OR_RETURN(total_file_size, file_w->FileSize());

  CHECK(bytes_read <= total_file_size);
  if (total_file_size != bytes_read) {
    LOG(ERROR) << "Unable to consume entire chunk file, only read "
               << entries.size() << " entries";
  }

  KeySet keys;
  for (const std::string& key : header.keys()) {
    keys.insert(key);
  }

  LOG(INFO) << "Loaded chunk with " << entries.size() << " entries in "
            << timer.TimeSoFarMillis().count() << "ms";
  return make_unique<Chunk>(file_opener, chunk_id, keys, &entries);
}

StatusOr<std::unique_ptr<Chunk>> Chunk::InitChunk(FileOpener* file_opener,
                                                  const KeySet& keys,
                                                  uint64_t chunk_id) {
  FWrapper* file_w;
  ASSIGN_OR_RETURN(file_w, file_opener->GetStream(chunk_id));

  ChunkFileHeaderProto header;
  for (const std::string& key : keys) {
    header.add_keys(key);
  }

  uint64_t header_size = header.ByteSize();
  RETURN_IF_ERROR(file_w->WriteUint32(header_size));

  std::string serialized = header.SerializeAsString();
  RETURN_IF_ERROR(file_w->Write(serialized.data(), serialized.size()));

  return make_unique<Chunk>(file_opener, chunk_id, keys);
}

StatusOr<const EntryHandle*> Chunk::AppendEntry(
    const std::vector<Value>& values, const Entry& entry) {
  using namespace std::chrono;
  CHECK(values.size() == keys_.size()) << "values " << values.size() << " vs "
                                       << " keys " << keys_.size();

  FWrapper* chunk_file_w;
  ASSIGN_OR_RETURN(chunk_file_w, file_opener_->GetStream(chunk_id_));

  uint64_t file_size;
  ASSIGN_OR_RETURN(file_size, chunk_file_w->Tell());

  auto new_handle = make_unique<EntryHandle>(values, TimeNow(), chunk_id_,
                                             file_size, entry.size());
  EntryHandleProto proto = new_handle->ToProto();
  std::string serialized = proto.SerializeAsString();

  RETURN_IF_ERROR(chunk_file_w->WriteUint32(serialized.size()));
  RETURN_IF_ERROR(chunk_file_w->Write(serialized.data(), serialized.size()));
  RETURN_IF_ERROR(chunk_file_w->Write(&(entry[0]), entry.size()));

  const EntryHandle* ptr = new_handle.get();
  entries_.emplace_back(std::move(new_handle));
  return ptr;
}

StatusOr<Entry> Chunk::ReadEntry(const EntryHandle& handle) {
  uint64_t entry_size = handle.entry_size();

  FWrapper* chunk_file_w;
  ASSIGN_OR_RETURN(chunk_file_w, file_opener_->GetStream(chunk_id_));
  RETURN_IF_ERROR(chunk_file_w->Seek(handle.offset_in_chunk()));

  uint32_t handle_size;
  ASSIGN_OR_RETURN(handle_size, chunk_file_w->ReadUint32());
  chunk_file_w->Skip(handle_size);

  Entry out;
  out.resize(entry_size);
  RETURN_IF_ERROR(chunk_file_w->Read(entry_size, &(out[0])));

  return out;
}

}  // namespace deep_storage
}  // namespace nc
