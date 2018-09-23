#include "deep_storage.h"

#include <gtest/gtest.h>
#include <ncode/file.h>
#include <ncode/port.h>
#include <ncode/status.h>
#include <ncode/strutil.h>
#include <string>

namespace nc {
namespace deep_storage {
namespace test {

class FileOpenerFixture : public ::testing::Test {
 protected:
  void SetUp() override { File::DeleteRecursively("f_1", nullptr, nullptr); }
};

TEST_F(FileOpenerFixture, BadOpen) {
  ASSERT_DEATH(FileOpener(0, "f_$0"), ".*");
}

TEST_F(FileOpenerFixture, Open) {
  FileOpener opener(1, "f_$0");

  auto result = opener.GetStream(1);
  ASSERT_TRUE(result.ok());
  ASSERT_NE(nullptr, result.ValueOrDie());
}

TEST_F(FileOpenerFixture, OpenAppend) {
  FileOpener file_opener(1, "f_$0");

  FWrapper* file_w = file_opener.GetStream(1).ValueOrDie();

  ASSERT_TRUE(file_w->WriteUint32(42UL).ok());
  ASSERT_TRUE(file_w->Seek(0).ok());
  ASSERT_EQ(42UL, file_w->ReadUint32().ValueOrDie());
}

class ChunkFixture : public ::testing::Test {
 protected:
  ChunkFixture() : file_opener_(100, "f_$0") {}

  void SetUp() override {
    for (uint64_t i = 0; i < 100; ++i) {
      File::DeleteRecursively(Substitute("f_$0", i), nullptr, nullptr);
    }
  }

  FileOpener file_opener_;
};

TEST_F(ChunkFixture, InitAndRead) {
  KeySet keys;
  keys.insert("A");

  Chunk::InitChunk(&file_opener_, keys, 1).ConsumeValueOrDie();
  auto chunk_two = Chunk::ReadChunk(&file_opener_, 1).ConsumeValueOrDie();

  ASSERT_EQ(chunk_two->keys(), keys);
}

TEST_F(ChunkFixture, InitAndBadAppend) {
  KeySet keys;
  keys.insert("Integer key");

  auto chunk_one = Chunk::InitChunk(&file_opener_, keys, 1).ConsumeValueOrDie();

  Entry entry;
  std::vector<Value> values = {static_cast<uint64_t>(1), "second_value"};

  ASSERT_DEATH(chunk_one->AppendEntry(values, entry).ValueOrDie(), ".*");
}

TEST_F(ChunkFixture, InitAndAppend) {
  KeySet keys;
  keys.insert("Integer key");
  keys.insert("String key");

  auto chunk_one = Chunk::InitChunk(&file_opener_, keys, 1).ConsumeValueOrDie();
  Entry entry = {1u, 2u, 3u};
  std::vector<Value> values = {static_cast<uint64_t>(1), "second_value"};
  auto result = chunk_one->AppendEntry(values, entry);
  ASSERT_TRUE(result.ok());

  auto chunk_two = Chunk::InitChunk(&file_opener_, keys, 1).ConsumeValueOrDie();
  auto read_status = chunk_two->ReadEntry(*(result.ValueOrDie()));
  ASSERT_TRUE(read_status.ok());
  ASSERT_EQ(entry, read_status.ValueOrDie());
}

TEST_F(ChunkFixture, MultiAppend) {
  KeySet keys;
  keys.insert("Integer key");
  keys.insert("String key");

  std::vector<const EntryHandle*> handles;
  Entry entry = {1u, 2u, 3u};
  auto chunk_one = Chunk::InitChunk(&file_opener_, keys, 1).ConsumeValueOrDie();
  for (uint32_t i = 0; i < 1000000u; ++i) {
    std::vector<Value> values = {static_cast<uint64_t>(i), "second_value"};
    auto result = chunk_one->AppendEntry(values, entry);
    ASSERT_TRUE(result.ok());
    handles.emplace_back(result.ValueOrDie());
  }

  {
    auto chunk_two = Chunk::ReadChunk(&file_opener_, 1).ConsumeValueOrDie();
    for (uint32_t i = 0; i < 1000000u; ++i) {
      auto read_result = chunk_two->ReadEntry(*handles[i]);
      ASSERT_TRUE(read_result.ok());
      ASSERT_EQ(entry, read_result.ValueOrDie());
    }
  }

  using namespace std::chrono;
  auto stats = file_opener_.GetStream(1).ValueOrDie()->stats();
  LOG(INFO) << "read "
            << duration_cast<milliseconds>(stats.total_read_time).count()
            << "ms, write "
            << duration_cast<milliseconds>(stats.total_write_time).count()
            << "ms, seek "
            << duration_cast<milliseconds>(stats.total_seek_time).count()
            << "ms";
}

// TEST(Chunk, AppendOver) {
//  Chunk chunk_one = Chunk::InitChunk(kTestFile, 1000).ValueOrDie();
//
//  for (uint32_t i = 0; i < 1000u; ++i) {
//    std::string some_value = StrCat("some_value_", i);
//    auto status = chunk_one.AppendEntry(&(some_value[0]), some_value.size());
//    ASSERT_TRUE(status.ok());
//  }
//
//  std::string some_value = "some_value";
//  auto status = chunk_one.AppendEntry(&(some_value[0]), some_value.size());
//  ASSERT_FALSE(status.ok());
//}

}  // namespace test
}  // namespace deep_storage
}  // namespace nc
