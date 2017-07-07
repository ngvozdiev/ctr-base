#include <gflags/gflags.h>
#include <chrono>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ncode_common/src/logging.h"
#include "ncode_common/src/map_util.h"
#include "ncode_common/src/strutil.h"
#include "metrics/metrics_parser.h"

DEFINE_string(sp_input, "", "The sp metrics file.");
DEFINE_string(minmax_input, "", "The minmax metrics file.");
DEFINE_string(metric, "tcp_source_completion_times", "The metric.");

using namespace nc;
using namespace nc::metrics::parser;

using DataVector = std::vector<std::pair<uint64_t, double>>;

//void DumpDataVector(const DataVector& data_vector) {
//  std::vector<std::string> out;
//  for (const auto& ts_and_value : data_vector) {
//    out.emplace_back(
//        nc::StrCat("(", ts_and_value.first, ", ", ts_and_value.second, ")"));
//  }
//
//  out.resize(std::min(out.size(), 20ul));
//  LOG(ERROR) << nc::Join(out, ",");
//}

static std::vector<std::chrono::milliseconds> GetDeltas(
    const DataVector& sp_data, const DataVector& other_data) {
  size_t min_data_size = std::min(sp_data.size(), other_data.size());
  std::vector<std::chrono::milliseconds> out(min_data_size);
  for (size_t i = 0; i < min_data_size; ++i) {
    uint64_t delta = other_data[i].second - sp_data[i].second;
    CHECK(delta >= 0);
    out[i] = std::chrono::milliseconds(delta);
  }

  return out;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_sp_input.empty()) << "need --sp_input";
  CHECK(!FLAGS_minmax_input.empty()) << "need --sp_input";

  std::map<std::pair<std::string, std::string>, DataVector> sp_data =
      SimpleParseNumericData(FLAGS_sp_input, FLAGS_metric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);

  std::map<std::pair<std::string, std::string>, DataVector> minmax_data =
      SimpleParseNumericData(FLAGS_minmax_input, FLAGS_metric, ".*", 0,
                             std::numeric_limits<uint64_t>::max(), 0);
  CHECK(sp_data.size() == minmax_data.size());

  std::vector<std::chrono::milliseconds> all_values;
  for (auto& id_and_data : sp_data) {
    const std::pair<std::string, std::string>& id = id_and_data.first;
    DataVector& data = id_and_data.second;
    DataVector& other_data = nc::FindOrDieNoPrint(minmax_data, id);
    CHECK(id.first == FLAGS_metric);

    std::vector<std::chrono::milliseconds> values = GetDeltas(data, other_data);
    all_values.insert(all_values.end(), values.begin(), values.end());
  }
}
