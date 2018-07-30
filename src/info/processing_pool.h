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
bool BlockingWriteProtoMessage(const google::protobuf::Message& message,
                               int socket);

// Performs a blocking read of a message from a socket.
void BlockingReadProtoMessage(int socket, google::protobuf::Message* message);

// Gets a timestamp.
std::chrono::milliseconds TimeNow();

template <typename Input>
using Processor = std::function<void(std::unique_ptr<Input>, uint32_t)>;

template <typename Input>
using InputQueue = nc::PtrQueue<Input, 1024>;

template <typename Input>
class ProcessingPool {
 public:
  ProcessingPool(Processor<Input> processor, InputQueue<Input>* input_queue,
                 size_t thread_count)
      : processor_(processor), input_queue_(input_queue) {
    for (size_t i = 0; i < thread_count; ++i) {
      threads_.emplace_back([this, i] { RunWorker(i); });
    }
  }

  ~ProcessingPool() {
    to_kill_ = true;
    for (size_t i = 0; i < threads_.size(); ++i) {
      threads_[i].join();
    }
  }

 private:
  void RunWorker(uint32_t worker_index) {
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

      processor_(std::move(input), worker_index);
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
  ProtobufServer(const nc::viz::TCPServerConfig& server_config,
                 size_t thread_count,
                 std::chrono::milliseconds log_message_timeout =
                     std::chrono::milliseconds(10))
      : thread_count_(thread_count),
        log_message_timeout_(log_message_timeout),
        tcp_server_(server_config, ParseConstantSizeHeader, &input_queue_) {}

  virtual ~ProtobufServer() { Stop(); }

  void Start() {
    processing_pool_ =
        nc::make_unique<ProcessingPool<nc::viz::IncomingHeaderAndMessage>>(
            [this](std::unique_ptr<nc::viz::IncomingHeaderAndMessage>
                       header_and_message,
                   uint32_t worker_index) {
              HandleIncomingMessage(std::move(header_and_message),
                                    worker_index);
            },
            &input_queue_, thread_count_);
    log_processing_pool_ = nc::make_unique<ProcessingPool<LoggedMessage>>(
        [this](std::unique_ptr<LoggedMessage> logged_message,
               uint32_t worker_index) {
          nc::Unused(worker_index);
          LogMessage(*logged_message);
        },
        &log_queue_, 1);

    tcp_server_.Start();
  }

  void Stop() { tcp_server_.Stop(); }

  void Join() { tcp_server_.Join(); }

 protected:
  struct LoggedMessage {
    LoggedMessage(
        std::unique_ptr<nc::viz::IncomingHeaderAndMessage> header_and_message,
        std::unique_ptr<Request> request, uint32_t worker_index,
        uint32_t response_size_bytes,
        std::chrono::milliseconds processing_started,
        std::chrono::milliseconds processing_done,
        std::chrono::milliseconds response_sent)
        : header_and_message(std::move(header_and_message)),
          request(std::move(request)),
          worker_index(worker_index),
          response_size_bytes(response_size_bytes),
          processing_started(processing_started),
          processing_done(processing_done),
          response_sent(response_sent) {}

    std::unique_ptr<nc::viz::IncomingHeaderAndMessage> header_and_message;
    std::unique_ptr<Request> request;
    uint32_t worker_index;
    uint32_t response_size_bytes;

    std::chrono::milliseconds processing_started;
    std::chrono::milliseconds processing_done;
    std::chrono::milliseconds response_sent;
  };

  virtual std::unique_ptr<Reply> HandleRequest(const Request& request) = 0;

  virtual void LogMessage(const LoggedMessage& logged_message) {
    nc::Unused(logged_message);
  }

 private:
  using LogQueue = PtrQueue<LoggedMessage, 1024>;

  void HandleIncomingMessage(
      std::unique_ptr<nc::viz::IncomingHeaderAndMessage> header_and_message,
      uint32_t worker_index) {
    size_t message_size =
        header_and_message->buffer.size() - header_and_message->header_offset;
    auto request = nc::make_unique<Request>();
    request->ParseFromArray(
        header_and_message->buffer.data() + header_and_message->header_offset,
        message_size);

    std::chrono::milliseconds processing_started = TimeNow();
    std::unique_ptr<Reply> reply = HandleRequest(*request);
    std::chrono::milliseconds processing_done = TimeNow();

    uint32_t response_size_bytes = reply->ByteSize();
    if (reply) {
      BlockingWriteProtoMessage(*reply, header_and_message->socket);
    }
    std::chrono::milliseconds response_sent = TimeNow();

    auto logged_message = nc::make_unique<LoggedMessage>(
        std::move(header_and_message), std::move(request), worker_index,
        response_size_bytes, processing_started, processing_done,
        response_sent);

    bool timed_out;
    log_queue_.ProduceOrBlockWithTimeout(std::move(logged_message),
                                         log_message_timeout_, &timed_out);
    if (timed_out) {
      LOG(INFO) << "Logging too slow; unable to log message";
    }
  }

  const size_t thread_count_;
  const std::chrono::milliseconds log_message_timeout_;
  nc::viz::IncomingMessageQueue input_queue_;
  LogQueue log_queue_;
  std::unique_ptr<ProcessingPool<nc::viz::IncomingHeaderAndMessage>>
      processing_pool_;
  std::unique_ptr<ProcessingPool<LoggedMessage>> log_processing_pool_;
  nc::viz::TCPServer tcp_server_;
};

}  // namespace nc
#endif
