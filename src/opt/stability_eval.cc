#include <gflags/gflags.h>
#include <stddef.h>
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/file.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/lp/demand_matrix.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/net/net_gen.h"
#include "ncode_common/src/strutil.h"
#include "../common.h"
#include "ctr.h"
#include "path_provider.h"

DEFINE_string(topology_file, "", "A topology file");
DEFINE_string(matrices, "", "One or more matrices");
DEFINE_double(demand_fraction, 0.05, "By how much to vary demand");
DEFINE_double(flow_count_fraction, std::numeric_limits<double>::max(),
              "By how much to vary flow counts. If set to max will update flow "
              "counts to be proportional to demand");
DEFINE_uint64(aggregate_count, 1, "How many aggregates to change");
DEFINE_uint64(try_count, 1, "How many tries to perform");
DEFINE_string(opt, "CTR", "Optimizer to run. One of CTR,MinMax,B4");
DEFINE_uint64(seed, 1ul, "Seed to use for the RNG");

// Each field has a datatype and an identifier.
using HeaderField = std::pair<std::string, std::string>;

template <typename... Args>
struct HeaderWriter {
  static void Execute(std::vector<std::string>* out) { nc::Unused(out); }
};

template <typename T>
struct HeaderWriter<T> {
  static void Execute(std::vector<std::string>* out) {
    LOG(FATAL) << "Unknown type";
    nc::Unused(out);
  }
};

template <>
struct HeaderWriter<std::string> {
  static void Execute(std::vector<std::string>* out) {
    out->emplace_back("S256");
  }
};

template <>
struct HeaderWriter<uint8_t> {
  static void Execute(std::vector<std::string>* out) {
    out->emplace_back("u1");
  }
};

template <>
struct HeaderWriter<uint16_t> {
  static void Execute(std::vector<std::string>* out) {
    out->emplace_back("u2");
  }
};

template <>
struct HeaderWriter<uint32_t> {
  static void Execute(std::vector<std::string>* out) {
    out->emplace_back("u4");
  }
};

template <>
struct HeaderWriter<uint64_t> {
  static void Execute(std::vector<std::string>* out) {
    out->emplace_back("u8");
  }
};

template <>
struct HeaderWriter<float> {
  static void Execute(std::vector<std::string>* out) {
    out->emplace_back("f4");
  }
};

template <>
struct HeaderWriter<double> {
  static void Execute(std::vector<std::string>* out) {
    out->emplace_back("f8");
  }
};

template <typename T, typename... Args>
struct HeaderWriter<T, Args...> {
  static void Execute(std::vector<std::string>* out) {
    HeaderWriter<T> current;
    current.Execute(out);

    HeaderWriter<Args...> rest;
    rest.Execute(out);
  }
};

// static void WriteU16(uint16_t value, std::vector<uint8_t>* out) {
//  out->emplace_back(static_cast<uint8_t>(value));
//  out->emplace_back(static_cast<uint8_t>(value >> 8));
//}
//
// static void WriteU32(uint16_t value, std::vector<uint8_t>* out) {
//  out->emplace_back(static_cast<uint8_t>(value));
//  out->emplace_back(static_cast<uint8_t>(value >> 8));
//  out->emplace_back(static_cast<uint8_t>(value >> 16));
//  out->emplace_back(static_cast<uint8_t>(value >> 24));
//}
//
// static void WriteU64(uint16_t value, std::vector<uint8_t>* out) {
//  out->emplace_back(static_cast<uint8_t>(value));
//  out->emplace_back(static_cast<uint8_t>(value >> 8));
//  out->emplace_back(static_cast<uint8_t>(value >> 16));
//  out->emplace_back(static_cast<uint8_t>(value >> 24));
//  out->emplace_back(static_cast<uint8_t>(value >> 32));
//  out->emplace_back(static_cast<uint8_t>(value >> 40));
//  out->emplace_back(static_cast<uint8_t>(value >> 48));
//  out->emplace_back(static_cast<uint8_t>(value >> 56));
//}
//

template <typename... Args>
class NpyWriter {
 public:
  explicit NpyWriter(const std::vector<std::string>& names) {
    std::vector<std::string> data_types;
    HeaderWriter<Args...>::Execute(&data_types);

    CHECK(data_types.size() == names.size());
    std::vector<std::string> dtypes_as_str;
    for (size_t i = 0; i < data_types.size(); ++i) {
      dtypes_as_str.emplace_back(
          nc::StrCat("('", names[i], "','", data_types[i], "')"));
    }

    dtype_ = nc::StrCat("[", nc::Join(dtypes_as_str, ","), "]");
  }

  const std::string& dtype() const { return dtype_; }

  template <class... Vs>
  void AddData(Vs&&... vals) {
    data_.emplace_back(vals...);
  }

  std::string Serialize() const {
    std::string dict = nc::StrCat("{'descr' : \"", dtype_,
                                  "\", 'fortran_order': False, 'shape': (",
                                  std::to_string(data_.size()), ",)} ");
    int remainder = 16 - (10 + dict.size()) % 16;
    dict.insert(dict.end(), remainder, ' ');
    dict.back() = '\n';

    std::string header = nc::StrCat(
        static_cast<uint8_t>(0x93), "NUMPY", static_cast<uint8_t>(0x01),
        static_cast<uint8_t>(0x00), static_cast<uint16_t>(dict.size()));
    header.insert(header.end(), dict.begin(), dict.end());
    return header;
  }

 private:
  std::string dtype_;
  std::vector<std::tuple<Args...>> data_;
};

template <typename T>
static std::string PercentilesToString(
    std::vector<T>* values,
    const std::vector<size_t>& percenile_indices = {50, 90, 95, 100}) {
  std::vector<T> percentiles = nc::Percentiles(values);
  std::vector<std::string> out;
  for (size_t p : percenile_indices) {
    CHECK(p <= 100);
    out.emplace_back(std::to_string(percentiles[p]));
  }

  return nc::Join(out, " ");
}

static void ParseMatrix(const std::string& tm_id, const ctr::TrafficMatrix& tm,
                        size_t seed, ctr::Optimizer* optimizer) {
  std::unique_ptr<ctr::RoutingConfiguration> baseline = optimizer->Optimize(tm);
  std::vector<ctr::AggregateDelta> deltas;
  double worst_total_demand_delta = 0;
  double worst_total_flow_delta = 0;
  //  size_t worst_demand_i = 0;
  std::vector<double> fraction_deltas;
  std::vector<size_t> route_add_counts;
  std::vector<size_t> route_update_counts;
  std::vector<size_t> route_remove_counts;

  //  LOG(INFO) << baseline->ToString();

  for (size_t i = 0; i < FLAGS_try_count; ++i) {
    std::mt19937 rnd(seed + i);
    auto rnd_tm = tm.Randomize(FLAGS_demand_fraction, FLAGS_flow_count_fraction,
                               FLAGS_aggregate_count, &rnd);
    if (!rnd_tm->ToDemandMatrix()->IsFeasible({})) {
      continue;
    }

    auto output = optimizer->Optimize(*rnd_tm);
    //    LOG(INFO) << output->ToString();
    ctr::RoutingConfigurationDelta routing_config_delta =
        baseline->GetDifference(*output);
    //    for (const auto& aggregate_and_delta :
    //    routing_config_delta.aggregates) {
    //      if (aggregate_and_delta.second.fraction_delta > 0.01) {
    //        LOG(INFO) << "C " <<
    //        aggregate_and_delta.first.ToString(*tm.graph())
    //                  << " " << aggregate_and_delta.second.fraction_delta << "
    //                  "
    //                  << aggregate_and_delta.second.path_stretch_gain;
    //      }
    //    }

    double demand_delta = routing_config_delta.total_volume_fraction_delta;
    double flow_delta = routing_config_delta.total_flow_fraction_delta;

    //    if (demand_delta > worst_total_demand_delta) {
    //      worst_demand_i = i;
    //    }

    worst_total_demand_delta = std::max(worst_total_demand_delta, demand_delta);
    worst_total_flow_delta = std::max(worst_total_flow_delta, flow_delta);

    for (const auto& aggregate_id_and_delta : routing_config_delta.aggregates) {
      const ctr::AggregateDelta& aggregate_delta =
          aggregate_id_and_delta.second;
      fraction_deltas.emplace_back(aggregate_delta.fraction_delta);
    }

    size_t route_adds;
    size_t route_removals;
    size_t route_updates;
    std::tie(route_adds, route_removals, route_updates) =
        routing_config_delta.TotalRoutes();

    route_add_counts.emplace_back(route_adds);
    route_remove_counts.emplace_back(route_removals);
    route_update_counts.emplace_back(route_updates);
  }

  std::cout << tm_id << " " << worst_total_demand_delta << " "
            << tm.demands().size() << " " << worst_total_flow_delta << " "
            << PercentilesToString(&fraction_deltas) << " "
            << PercentilesToString(&route_add_counts) << " "
            << PercentilesToString(&route_update_counts) << " "
            << PercentilesToString(&route_remove_counts) << "\n";
}

// The TM that we load will have no flow counts. Need some out of thin air.
static std::map<nc::lp::SrcAndDst, size_t> GetFlowCountMap(
    const nc::lp::DemandMatrix& demand_matrix) {
  std::map<nc::lp::SrcAndDst, size_t> out;
  for (const auto& element : demand_matrix.elements()) {
    out[{element.src, element.dst}] = 1000;
  }

  return out;
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_topology_file.empty());

  std::vector<std::string> matrix_files;
  std::vector<std::string> to_glob = nc::Split(FLAGS_matrices, ",");
  for (const std::string string_to_glob : to_glob) {
    std::vector<std::string> globbed = nc::Glob(string_to_glob);
    matrix_files.insert(matrix_files.end(), globbed.begin(), globbed.end());
  }
  CHECK(!matrix_files.empty());

  std::vector<std::string> nodes_in_order;
  nc::net::GraphBuilder graph_builder = nc::net::LoadRepetitaOrDie(
      nc::File::ReadFileToStringOrDie(FLAGS_topology_file), &nodes_in_order);
  graph_builder.FuzzLinkDelays(std::chrono::microseconds(100));
  nc::net::GraphStorage graph(graph_builder);

  for (const std::string& matrix_file : matrix_files) {
    auto demand_matrix = nc::lp::DemandMatrix::LoadRepetitaOrDie(
        nc::File::ReadFileToStringOrDie(matrix_file), nodes_in_order, &graph);
    ctr::TrafficMatrix traffic_matrix(*demand_matrix,
                                      GetFlowCountMap(*demand_matrix));

    auto path_provider = nc::make_unique<ctr::PathProvider>(&graph);
    std::unique_ptr<ctr::Optimizer> optimizer;
    if (FLAGS_opt == "CTR") {
      optimizer = nc::make_unique<ctr::CTROptimizer>(std::move(path_provider));
    } else if (FLAGS_opt == "MinMax") {
      optimizer =
          nc::make_unique<ctr::MinMaxOptimizer>(std::move(path_provider));
    } else if (FLAGS_opt == "B4") {
      optimizer = nc::make_unique<ctr::B4Optimizer>(std::move(path_provider));
    } else {
      LOG(FATAL) << "Unknown optimizer";
    }

    ParseMatrix(nc::StrCat(FLAGS_topology_file, ":", matrix_file),
                traffic_matrix, FLAGS_seed, optimizer.get());
  }
}
