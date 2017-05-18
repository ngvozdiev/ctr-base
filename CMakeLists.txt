cmake_minimum_required(VERSION 2.8.7)

# If there is already a target called ncode_common someone has already
# added this project. Will ignore.
if(TARGET ctr_base)
  return()
endif(TARGET ctr_base)

project(NCodeLP)
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

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_extensions)

if (NOT TARGET ncode_common)
   EXECUTE_PROCESS(COMMAND "git" "clone" "https://github.com/ngvozdiev/ncode-common.git" "${CMAKE_BINARY_DIR}/external/ncode_common")
   add_subdirectory(${CMAKE_BINARY_DIR}/external/ncode_common)
endif()

if (NOT TARGET ncode_net)
   EXECUTE_PROCESS(COMMAND "git" "clone" "https://github.com/ngvozdiev/ncode-net.git" "${CMAKE_BINARY_DIR}/external/ncode_net")
   add_subdirectory(${CMAKE_BINARY_DIR}/external/ncode_net)
endif()

if (NOT TARGET ncode_lp)
   EXECUTE_PROCESS(COMMAND "git" "clone" "https://github.com/ngvozdiev/ncode-lp.git" "${CMAKE_BINARY_DIR}/external/ncode_lp")
   add_subdirectory(${CMAKE_BINARY_DIR}/external/ncode_lp)
endif()

include_directories(${CMAKE_SOURCE_DIR}/external ${PROJECT_BINARY_DIR} ${OPTIMIZER_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/external)

set(CTR_BASE_HEADER_FILES src/opt/opt.h src/opt/path_provider.h src/common.h src/opt/ctr.h)
add_library(ctr_base STATIC src/opt/opt.cc src/opt/path_provider.cc src/common.cc src/opt/ctr.cc ${CTR_BASE_HEADER_FILES})
target_link_libraries(ctr_base ncode_lp ${OPTIMIZER_LIBRARIES})

if (NOT CTR_BASE_DISABLE_TESTS)
   add_test_exec(ctr_opt_test src/opt/opt_test.cc ctr_base)
endif()
