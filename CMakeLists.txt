cmake_minimum_required(VERSION 2.8.7)

# If there is already a target called ncode_base someone has already
# added this project. Will ignore.
if(TARGET ctr_base)
  return()
endif(TARGET ctr_base)

project(CTR)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_extensions)
option(CTR_BASE_DISABLE_TESTS "If tests should be compiled or not" ON)
option(CTR_BASE_DEBUG "A debug build" OFF)
option(CTR_BASE_ASAN "Compile with ASAN on" OFF)
option(CTR_BASE_TSAN "Compile with TSAN on" OFF)

set(CTR_BASE_BASE_FLAGS "-g -std=c++11 -pedantic-errors -Winit-self -Woverloaded-virtual -Wuninitialized -Wall -Wextra -fno-exceptions")
set(CTR_BASE_BASE_LD_FLAGS "")
if (CTR_BASE_ASAN)
   set(CTR_BASE_BASE_FLAGS "${CTR_BASE_BASE_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")
   set(CTR_BASE_BASE_LD_FLAGS "${CTR_BASE_BASE_LD_FLAGS} -fsanitize=address")
endif()
if (CTR_BASE_TSAN)
   set(CTR_BASE_BASE_FLAGS "${CTR_BASE_BASE_FLAGS} -fsanitize=thread -fno-omit-frame-pointer -fno-optimize-sibling-calls")
   set(CTR_BASE_BASE_LD_FLAGS "${CTR_BASE_BASE_LD_FLAGS} -fsanitize=thread")
endif()

if (CTR_BASE_DEBUG)
  set(CTR_BASE_BASE_FLAGS "${CTR_BASE_BASE_FLAGS} -O0 -fno-omit-frame-pointer")
else()
  set(CTR_BASE_BASE_FLAGS "${CTR_BASE_BASE_FLAGS} -O3 -march=native -DNDEBUG")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CTR_BASE_BASE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${CTR_BASE_BASE_LD_FLAGS}")

if (NOT CTR_BASE_DISABLE_TESTS)
   include(CTest)
   add_subdirectory(external/googletest)
   macro(add_test_exec name src_file deps)
     add_executable(${name} ${src_file})
     target_link_libraries(${name} gtest gmock_main ${deps} ${ARGN})
     add_test(NAME ${name} COMMAND ${name})
   endmacro(add_test_exec)
endif()

SET(OLD_BUILD_TESTING ${BUILD_TESTING})
SET(BUILD_TESTING False)
add_subdirectory(external/gflags)
SET(BUILD_TESTING ${OLD_BUILD_TESTING})

if (NOT TARGET ncode_common)
   EXECUTE_PROCESS(COMMAND "git" "clone" "https://github.com/ngvozdiev/ncode-common.git" "${CMAKE_BINARY_DIR}/external/ncode_common")
   add_subdirectory(${CMAKE_BINARY_DIR}/external/ncode_common)
endif()

find_package(Protobuf REQUIRED)
find_package(FFTW REQUIRED)
include_directories(${PROTOBUF_INCLUDE_DIRS} ${CMAKE_BINARY_DIR} ${FFTW_INCLUDES})
include_directories(${CMAKE_SOURCE_DIR}/external ${CMAKE_BINARY_DIR}/external)

PROTOBUF_GENERATE_CPP(PROTO_METRICS_SRCS PROTO_METRICS_HDRS src/metrics/metrics.proto)
set_property(SOURCE ${PROTO_METRICS_SRCS} APPEND_STRING PROPERTY COMPILE_FLAGS "-Wno-extended-offsetof")
set(METRICS_HEADER_FILES src/metrics/metrics.h src/metrics/metrics_parser.h ${PROTO_METRICS_HDRS})
add_library(metrics STATIC src/metrics/metrics.cc ${PROTO_METRICS_SRCS} ${METRICS_HEADER_FILES})
target_link_libraries(metrics ${PROTOBUF_LIBRARIES} ncode gflags)

add_library(metrics_parser src/metrics/metrics_parser.cc ${PROTO_METRICS_SRCS} ${METRICS_HEADER_FILES})
target_link_libraries(metrics_parser ${PROTOBUF_LIBRARIES} ncode)

add_library(geo_lib src/geo/geo.cc src/geo/geo.h)
target_link_libraries(geo_lib ncode)

PROTOBUF_GENERATE_CPP(PCAP_DATA_PROTO_SRCS PCAP_DATA_PROTO_HDRS src/pcap_data.proto)
set(CTR_BASE_HEADER_FILES src/opt/opt.h src/opt/path_provider.h src/common.h src/opt/ctr.h src/opt/oversubscription_model.h src/pcap_data.h ${PCAP_DATA_PROTO_HDRS} src/prob_model/dist_model.h src/prob_model/multi_rate_queue.h src/mean_est/mean_est.h src/routing_system.h src/net_mock.h src/tldr.h src/controller.h src/net_instrument.h src/opt_eval.h)
add_library(ctr_base STATIC src/opt/opt.cc src/opt/path_provider.cc src/common.cc src/opt/ctr.cc src/opt/oversubscription_model.cc src/pcap_data.cc ${PCAP_DATA_PROTO_SRCS} ${CTR_BASE_HEADER_FILES} src/prob_model/dist_model.cc src/prob_model/multi_rate_queue.cc src/mean_est/mean_est.cc src/routing_system.cc src/net_mock.cc src/controller.cc src/tldr.cc src/net_instrument.cc src/opt_eval.cc)
target_link_libraries(ctr_base ncode gflags_nothreads metrics ${FFTW_LIBRARIES} ${PROTOBUF_LIBRARIES})

add_executable(stability_eval src/opt/stability_eval.cc)
target_link_libraries(stability_eval ctr_base)

add_executable(stability_eval_dependent src/opt/stability_eval_dependent.cc)
target_link_libraries(stability_eval_dependent ctr_base)

add_executable(multi_rate_eval src/prob_model/multi_rate_eval.cc)
target_link_libraries(multi_rate_eval ctr_base)

add_executable(mean_est_eval src/mean_est/mean_est_eval.cc)
target_link_libraries(mean_est_eval ctr_base)

add_executable(pcap_data_util src/pcap_data_util)
target_link_libraries(pcap_data_util ctr_base)

add_executable(plot_pcap_data src/plot_pcap_data.cc)
target_link_libraries(plot_pcap_data ctr_base)

add_executable(net_mock_util src/net_mock_util)
target_link_libraries(net_mock_util ctr_base)

add_executable(ladder_tcp_util src/ladder_tcp_util)
target_link_libraries(ladder_tcp_util ctr_base)

add_executable(headroom_delay_util src/headroom_delay_util)
target_link_libraries(headroom_delay_util ctr_base)

if (NOT CTR_BASE_DISABLE_TESTS)
   add_test_exec(ctr_opt_test src/opt/opt_test.cc ctr_base)
   add_test_exec(dist_model_test src/prob_model/dist_model_test.cc ctr_base)

   add_library(metrics_test_util STATIC src/metrics/metrics_test_util.cc)
   target_link_libraries(metrics_test_util metrics gtest)

   add_test_exec(metrics_test src/metrics/metrics_test.cc metrics metrics_parser metrics_test_util)
   add_test_exec(metrics_parser_test src/metrics/metrics_parser_test.cc metrics metrics_parser metrics_test_util)

   add_test_exec(ctr_common_test src/common_test.cc ctr_base)
   
   add_test_exec(geo_test src/geo/geo_test.cc geo_lib)
endif()

add_executable(metrics_explore src/metrics/metrics_explore.cc)
target_link_libraries(metrics_explore metrics_parser gflags)

add_executable(tm_util src/tm_util.cc)
target_link_libraries(tm_util ctr_base gflags)

add_executable(meshy_util src/meshy_util.cc)
target_link_libraries(meshy_util ctr_base gflags)

add_executable(trace_pick_util src/trace_pick_util.cc)
target_link_libraries(trace_pick_util ctr_base gflags)

add_executable(plot_completion_time_delta src/plot_completion_time_delta.cc)
target_link_libraries(plot_completion_time_delta metrics_parser gflags)

add_executable(plot_queue_size_dist src/plot_queue_size_dist.cc)
target_link_libraries(plot_queue_size_dist metrics_parser ctr_base gflags)

add_executable(plot_algorithm_eval src/plot_algorithm_eval.cc)
target_link_libraries(plot_algorithm_eval metrics_parser gflags)

add_executable(plot_stability_eval src/plot_stability_eval.cc)
target_link_libraries(plot_stability_eval metrics_parser gflags)

add_executable(plot_top_grow src/plot_top_grow.cc)
target_link_libraries(plot_top_grow metrics_parser gflags)

add_executable(plot_headroom_delay src/plot_headroom_delay.cc)
target_link_libraries(plot_headroom_delay metrics_parser gflags)

add_executable(opt_eval_util src/opt_eval_util.cc)
target_link_libraries(opt_eval_util ctr_base)

add_library(flow_track_lib STATIC src/flow_track.cc)
target_link_libraries(flow_track_lib ncode gflags_nothreads)

add_executable(flow_track src/flow_track_main.cc)
target_link_libraries(flow_track flow_track_lib)

add_executable(geo_main src/geo/geo_main.cc)
target_link_libraries(geo_main geo_lib gflags)

add_executable(top_grow src/top_grow.cc)
target_link_libraries(top_grow geo_lib ctr_base metrics gflags)
