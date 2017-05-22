cmake_minimum_required(VERSION 2.8.7)

# If there is already a target called ncode_common someone has already
# added this project. Will ignore.
if(TARGET ctr_base)
  return()
endif(TARGET ctr_base)

project(CTR)
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

if (NOT TARGET ncode_common)
   EXECUTE_PROCESS(COMMAND "git" "clone" "https://github.com/ngvozdiev/ncode-common.git" "${CMAKE_BINARY_DIR}/external/ncode_common")
   add_subdirectory(${CMAKE_BINARY_DIR}/external/ncode_common)
endif()

include_directories(${CMAKE_SOURCE_DIR}/external ${CMAKE_BINARY_DIR}/external)
SET(OLD_BUILD_TESTING ${BUILD_TESTING})
SET(BUILD_TESTING False)
add_subdirectory(external/gflags)
SET(BUILD_TESTING ${OLD_BUILD_TESTING})

set(CTR_BASE_HEADER_FILES src/opt/opt.h src/opt/path_provider.h src/common.h src/opt/ctr.h src/opt/oversubscription_model.h)
add_library(ctr_base STATIC src/opt/opt.cc src/opt/path_provider.cc src/common.cc src/opt/ctr.cc src/opt/oversubscription_model.cc ${CTR_BASE_HEADER_FILES})
target_link_libraries(ctr_base ncode)

add_executable(stability_eval src/opt/stability_eval.cc)
target_link_libraries(stability_eval ctr_base gflags)

if (NOT CTR_BASE_DISABLE_TESTS)
   add_test_exec(ctr_opt_test src/opt/opt_test.cc ctr_base)
endif()
