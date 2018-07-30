#ifndef CTR_INFO_SERVER_H
#define CTR_INFO_SERVER_H

#include <memory>

#include "../../build/info.pb.h"
#include "info.h"
#include "processing_pool.h"

namespace ctr {

// Serves the information contained in InfoStorage.
class InfoServer
    : public nc::ProtobufServer<ctr::info::Request, ctr::info::Response> {
 public:
  InfoServer(std::unique_ptr<InfoStorage> info_storage);

  std::unique_ptr<ctr::info::Response> HandleRequest(
      const ctr::info::Request& request) override;

  void LogMessage(const LoggedMessage& logged_message) override;

 private:
  std::unique_ptr<ctr::info::Response> HandleSelect(
      const info::SelectInfoRequest& request) const;

  std::unique_ptr<ctr::info::Response> HandleTMGenRequest(
      const info::TrafficMatrixGenerateRequest& request) const;

  std::unique_ptr<InfoStorage> info_storage_;
};

}  // namespace ctr
#endif
