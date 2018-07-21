#ifndef CTR_INFO_SERVER_H
#define CTR_INFO_SERVER_H

#include <memory>

#include "../../build/info.pb.h"
#include "info.h"
#include "processing_pool.h"

namespace ctr {

// Serves the information contained in InfoStorage.
class InfoServer : public nc::ProtobufServer<ctr::info::InfoRequest,
                                             ctr::info::InfoResponse> {
 public:
  InfoServer(std::unique_ptr<InfoStorage> info_storage);

  std::unique_ptr<ctr::info::InfoResponse> HandleRequest(
      const ctr::info::InfoRequest& request) override;

 private:
  std::unique_ptr<ctr::info::InfoResponse> HandleTopologyInfo(
      const info::TopologyInfoRequest& request) const;

  std::unique_ptr<ctr::info::InfoResponse> HandleTrafficMatrixInfo(
      const info::TrafficMatrixInfoRequest& request) const;

  std::unique_ptr<ctr::info::InfoResponse> HandleRoutingInfo(
      const info::RoutingInfoRequest& request) const;

  std::unique_ptr<InfoStorage> info_storage_;
};

}  // namespace ctr
#endif
