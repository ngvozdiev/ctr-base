#include <gflags/gflags.h>
#include <cstdio>
#include <string>
#include <vector>

#include "ncode_common/src/logging.h"
#include "pcap_data.pb.h"
#include "pcap_data.h"

DEFINE_string(pcap_files, "", "Glob-string of .pcap files");
DEFINE_string(trace_location, "", "The location of the trace");
DEFINE_string(trace_date, "", "The date of the trace");
DEFINE_string(pcap_trace_store, "trace_store.pb",
              "A file with information about .pcap traces");
DEFINE_uint64(bin_size_microseconds, 1000, "Bin sizes in microseconds");
DEFINE_uint64(split_count, 100,
              "How many sub-traces to produce from the trace");
DEFINE_bool(summarize, true, "If true will print a summary of the store");

static void PrintSummary() {
  CHECK(!FLAGS_pcap_trace_store.empty());
  ctr::PcapTraceStore trace_store(FLAGS_pcap_trace_store);
  std::cout << trace_store.Summary();
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  if (FLAGS_summarize) {
    PrintSummary();
    return 0;
  }

  CHECK(!FLAGS_pcap_files.empty());
  CHECK(!FLAGS_trace_location.empty());
  CHECK(!FLAGS_trace_date.empty());

  std::vector<std::string> files = nc::Glob(FLAGS_pcap_files);
  CHECK(!files.empty());

  int day, month, year;
  const std::string& date = FLAGS_trace_date;
  CHECK(sscanf(date.c_str(), "%2d/%2d/%4d", &month, &day, &year) == 3);

  ctr::PBPcapDataTrace data_trace_pb;
  data_trace_pb.set_bin_size_micros(FLAGS_bin_size_microseconds);
  data_trace_pb.set_split_count(FLAGS_split_count);

  ctr::PBTraceId* trace_id = data_trace_pb.mutable_id();
  trace_id->set_location(FLAGS_trace_location);
  trace_id->set_day(day);
  trace_id->set_month(month);
  trace_id->set_year(year);

  for (const std::string& file : files) {
    data_trace_pb.add_pcap_files(file);
  }

  ctr::PcapDataTrace::Init(data_trace_pb, FLAGS_pcap_trace_store);
}
