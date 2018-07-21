#ifndef NCODE_PPOOL_H
#define NCODE_PPOOL_H

#include <ncode/ptr_queue.h>
#include <ncode/viz/server.h>
#include <stddef.h>
#include <atomic>
#include <chrono>
#include <functional>
#include <memory>
#include <thread>
#include <vector>

namespace google {
namespace protobuf {
class Message;
} /* namespace protobuf */
} /* namespace google */

namespace nc {

// A header parse function for fixed-size uint32_t headers.
bool ParseConstantSizeHeader(std::vector<char>::const_iterator from,
                             std::vector<char>::const_iterator to,
                             size_t* header_size, size_t* message_size);

// Performs a blocking send of a protobuf on a socket.
void BlockingWriteProtoMessage(const google::protobuf::Message& message,
                               int socket);

// Performs a blocking read of a message from a socket.
void BlockingReadProtoMessage(int socket, google::protobuf::Message* message);

template <typename Input>
using Processor = std::function<void(const Input&)>;

template <typename Input>
using InputQueue = nc::PtrQueue<Input, 1024>;

template <typename Input>
class ProcessingPool {
 public:
  ProcessingPool(Processor<Input> processor, InputQueue<Input>* input_queue,
                 size_t thread_count)
      : processor_(processor), input_queue_(input_queue) {
    for (size_t i = 0; i < thread_count; ++i) {
      threads_.emplace_back([this] { RunWorker(); });
    }
  }

  ~ProcessingPool() {
    to_kill_ = true;
    for (size_t i = 0; i < threads_.size(); ++i) {
      threads_[i].join();
    }
  }

 private:
  void RunWorker() {
    while (!to_kill_) {
      bool timed_out;
      std::unique_ptr<Input> input = input_queue_->ConsumeOrBlockWithTimeout(
          std::chrono::milliseconds(100), &timed_out);
      if (to_kill_) {
        break;
      }

      if (timed_out) {
        continue;
      }

      if (!input) {
        break;
      }

      processor_(*input);
    }
  }

  const Processor<Input> processor_;

  // Inputs come from here.
  InputQueue<Input>* input_queue_;

  // Whether or not to terminate the workers.
  std::atomic<bool> to_kill_;

  // The workers.
  std::vector<std::thread> threads_;
};

template <typename Request, typename Reply>
class ProtobufServer {
 public:
  ProtobufServer(uint32_t port, size_t thread_count)
      : thread_count_(thread_count),
        tcp_server_(port, ParseConstantSizeHeader, &input_queue_) {}

  virtual ~ProtobufServer() { Stop(); }

  virtual std::unique_ptr<Reply> HandleRequest(const Request& request) = 0;

  void Start() {
    processing_pool_ =
        nc::make_unique<ProcessingPool<nc::viz::IncomingHeaderAndMessage>>(
            [this](
                const nc::viz::IncomingHeaderAndMessage& header_and_message) {
              HandleIncomingMessage(header_and_message);
            },
            &input_queue_, thread_count_);
    tcp_server_.Start();
  }

  void Stop() { tcp_server_.Stop(); }

  void Join() { tcp_server_.Join(); }

 private:
  void HandleIncomingMessage(
      const nc::viz::IncomingHeaderAndMessage& header_and_message) {
    size_t message_size =
        header_and_message.buffer.size() - header_and_message.header_offset;
    Request request;
    request.ParseFromArray(
        header_and_message.buffer.data() + header_and_message.header_offset,
        message_size);

    std::unique_ptr<Reply> reply = HandleRequest(request);
    if (reply) {
      BlockingWriteProtoMessage(*reply, header_and_message.socket);
    }
  }

  const size_t thread_count_;
  nc::viz::IncomingMessageQueue input_queue_;
  std::unique_ptr<ProcessingPool<nc::viz::IncomingHeaderAndMessage>>
      processing_pool_;
  nc::viz::TCPServer tcp_server_;
};

}  // namespace nc
#endif
