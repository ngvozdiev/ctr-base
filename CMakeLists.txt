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
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${CTR_BASE_BASE_LD_FLAGS}")

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

PROTOBUF_GENERATE_CPP(PCAP_DATA_PROTO_SRCS PCAP_DATA_PROTO_HDRS src/pcap_data.proto)
set(CTR_BASE_HEADER_FILES src/opt/opt.h src/opt/path_provider.h src/common.h src/opt/ctr.h src/opt/oversubscription_model.h src/opt/opt_compare.h src/pcap_data.h src/flow_counter.h ${PCAP_DATA_PROTO_HDRS} src/prob_model/dist_model.h src/prob_model/multi_rate_queue.h src/mean_est/mean_est.h src/routing_system.h src/net_mock.h)
add_library(ctr_base STATIC src/opt/opt.cc src/opt/path_provider.cc src/common.cc src/opt/ctr.cc src/opt/oversubscription_model.cc src/opt/opt_compare.cc src/pcap_data.cc src/flow_counter.cc ${PCAP_DATA_PROTO_SRCS} ${CTR_BASE_HEADER_FILES} src/prob_model/dist_model.cc src/prob_model/multi_rate_queue.cc src/mean_est/mean_est.cc src/routing_system.cc src/net_mock.cc)
target_link_libraries(ctr_base ncode gflags_nothreads ${FFTW_LIBRARIES} ${PROTOBUF_LIBRARIES})

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

add_executable(net_mock_util src/net_mock_util)
target_link_libraries(net_mock_util ctr_base)

if (NOT CTR_BASE_DISABLE_TESTS)
   add_test_exec(ctr_opt_test src/opt/opt_test.cc ctr_base)
   add_test_exec(dist_model_test src/prob_model/dist_model_test.cc ctr_base)
endif()
