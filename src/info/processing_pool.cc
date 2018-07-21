#include "processing_pool.h"

#include <google/protobuf/message.h>
#include <ncode/logging.h>
#include <unistd.h>
#include <cstring>
#include <iterator>

namespace nc {

// A header parse function for fixed-size uint32_t headers.
bool ParseConstantSizeHeader(std::vector<char>::const_iterator from,
                             std::vector<char>::const_iterator to,
                             size_t* header_size, size_t* message_size) {
  if (std::distance(from, to) < 4) {
    return false;
  }

  *header_size = 4;
  uint32_t value = *(reinterpret_cast<const uint32_t*>(&(*from)));
  *message_size = ntohl(value);
  return true;
}

void BlockingWriteProtoMessage(const google::protobuf::Message& message,
                               int socket) {
  nc::viz::OutgoingHeaderAndMessage header_and_message(socket);

  std::vector<char>& buffer = header_and_message.buffer;
  uint32_t size = message.ByteSize();
  buffer.resize(size + 4);

  uint32_t header = htonl(size);
  memcpy(buffer.data(), &header, 4);
  CHECK(message.SerializeToArray(buffer.data() + 4, size))
      << "Unable to serialize message";
  nc::viz::BlockingWriteMessage(header_and_message);
}

void BlockingReadProtoMessage(int socket, google::protobuf::Message* message) {
  uint32_t header;
  CHECK(read(socket, &header, 4) == 4) << "Unable to read message header";

  uint32_t size = ntohl(header);
  std::vector<char> buffer(size);
  CHECK(read(socket, buffer.data(), size) == size) << "Unable to read message";
  CHECK(message->ParseFromArray(buffer.data(), size)) << "Invalid message";
}

}  // namespace nc
