#include <gflags/gflags.h>
#include <ncode/common.h>
#include <type_traits>

#include "info.h"
#include "info_server.h"
#include "processing_pool.h"

DEFINE_string(input, "info.pb", "All infos");

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  auto info_storage = nc::make_unique<ctr::InfoStorage>();
  info_storage->ReadFromFile(FLAGS_input);

  ctr::InfoServer info_server(std::move(info_storage));
  info_server.Start();
  info_server.Join();
}
