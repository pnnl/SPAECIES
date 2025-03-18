cmake_minimum_required(VERSION 3.22)

if( "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" )
  message( FATAL_ERROR
    " FATAL: In-source builds are not allowed."
    " Remember to delete CMakeCache.txt and any other generated files"
    " (e.g. the CMakeFiles/ directory), then build in a separate source tree.")
endif()

option(BUILD_TESTING "Enable testing." OFF)

project(SPAECIES LANGUAGES CXX VERSION 0.0.1)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
endif()

set(C_DBG_FLAGS "")
set(C_OPT_FLAGS -fno-math-errno -march=native)

set(CXX_DBG_FLAGS -Wall -fbounds-check)
set(CXX_OPT_FLAGS -fno-math-errno -march=native)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:${C_DBG_FLAGS}>")
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${CXX_DBG_FLAGS}>")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:${C_OPT_FLAGS}>")
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${CXX_OPT_FLAGS}>")
endif()

find_package(NetCDF REQUIRED)
find_package(SUNDIALS 7.2.1 QUIET COMPONENTS arkode sunlinsollapackdense)
if (NOT SUNDIALS_FOUND)
  include(FetchContent)
  set(FETCHCONTENT_QUIET OFF)
  FetchContent_Declare(sundials
    GIT_REPOSITORY https://github.com/LLNL/sundials.git
    GIT_TAG v7.2.1
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
  )
  set(BUILD_CVODE OFF)
  set(BUILD_CVODES OFF)
  set(BUILD_IDA OFF)
  set(BUILD_IDAS OFF)
  set(BUILD_KINSOL OFF)
  set(ENABLE_LAPACK ON)
  set(EXAMPLES_ENABLE_C OFF)
  set(BUILD_SHARED_LIBS OFF)
  FetchContent_MakeAvailable(sundials)
endif()
find_package(Boost REQUIRED COMPONENTS program_options)

if (BUILD_TESTING)
  find_package(Catch2 3 REQUIRED)
  include(CTest)
  include(Catch)
  add_subdirectory(tests)
endif()

add_subdirectory(libspaecies)
add_subdirectory(rainshaft)
