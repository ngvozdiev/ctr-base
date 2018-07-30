#include "processing_pool.h"

#include <gtest/gtest.h>
#include <ncode/common.h>

#include "../../build/info.pb.h"

namespace nc {
namespace test {

static constexpr viz::TCPServerConfig kConfig = {8080, 1 << 10};

class ProtobufServerForTest
    : public ProtobufServer<ctr::info::Request, ctr::info::Response> {
 public:
  ProtobufServerForTest() : ProtobufServer(kConfig, 10) {}

  std::unique_ptr<ctr::info::Response> HandleRequest(
      const ctr::info::Request& request) override {
    nc::Unused(request);
    return nc::make_unique<ctr::info::Response>();
  }
};

class ProtobufServerFixture : public ::testing::Test {
 protected:
  ProtobufServerForTest proto_server_;
};

TEST_F(ProtobufServerFixture, NoWork) {
  proto_server_.Start();
  proto_server_.Stop();
}

TEST_F(ProtobufServerFixture, NoWorkNoStop) { proto_server_.Stop(); }

TEST_F(ProtobufServerFixture, SingleMessage) {
  proto_server_.Start();
  int socket = nc::viz::Connect("127.0.0.1", 8080);

  ctr::info::Request request;
  BlockingWriteProtoMessage(request, socket);

  ctr::info::Response response;
  BlockingReadProtoMessage(socket, &response);
}

TEST_F(ProtobufServerFixture, MultiMessageSingleConnection) {
  proto_server_.Start();
  int socket = nc::viz::Connect("127.0.0.1", 8080);

  auto thread_one = std::thread([socket]() {
    ctr::info::Request request;
    for (size_t i = 0; i < 10000; ++i) {
      BlockingWriteProtoMessage(request, socket);
    }
  });

  auto thread_two = std::thread([socket]() {
    ctr::info::Response response;
    for (size_t i = 0; i < 10000; ++i) {
      BlockingReadProtoMessage(socket, &response);
    }
  });

  thread_one.join();
  thread_two.join();
}

TEST_F(ProtobufServerFixture, MultiConnectionMultiMessage) {
  proto_server_.Start();

  std::vector<int> sockets;
  for (size_t i = 0; i < 10; ++i) {
    int socket = nc::viz::Connect("127.0.0.1", 8080);
    sockets.emplace_back(socket);
  }

  auto thread_one = std::thread([&sockets]() {
    ctr::info::Request request;
    for (size_t i = 0; i < 10000; ++i) {
      int socket = sockets[i % 10];
      BlockingWriteProtoMessage(request, socket);
    }
  });

  auto thread_two = std::thread([&sockets]() {
    ctr::info::Response response;
    for (size_t i = 0; i < 10000; ++i) {
      int socket = sockets[i % 10];
      BlockingReadProtoMessage(socket, &response);
    }
  });

  thread_one.join();
  thread_two.join();
}

}  // namespace test
}  // namespace nc
