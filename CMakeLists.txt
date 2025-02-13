cmake_minimum_required(VERSION 3.22)

if( "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" )
  message( FATAL_ERROR
    " FATAL: In-source builds are not allowed."
    " Remember to delete CMakeCache.txt and any other generated files"
    " (e.g. the CMakeFiles/ directory), then build in a separate source tree.")
endif()

set(BUILD_TESTING OFF CACHE BOOL "Enable testing.")

project(SPAECIES LANGUAGES CXX VERSION 0.0.1)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(NetCDF REQUIRED)
find_package(SUNDIALS 7.2.0 QUIET COMPONENTS arkode sunlinsollapackdense)
if (NOT SUNDIALS_FOUND)
  include(FetchContent)
  set(FETCHCONTENT_QUIET OFF)
  FetchContent_Declare(sundials
    GIT_REPOSITORY https://github.com/LLNL/sundials.git
    GIT_TAG 0eff39663606f2ff280c4059a947ed62ae38180a # Release 7.2.0
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
  )
  set(BUILD_CVODE OFF)
  set(BUILD_CVODES OFF)
  set(BUILD_IDA OFF)
  set(BUILD_IDAS OFF)
  set(BUILD_KINSOL OFF)
  set(ENABLE_LAPACK ON)
  FetchContent_MakeAvailable(sundials)
endif()
find_package(Boost REQUIRED COMPONENTS program_options)

if (BUILD_TESTING)
  find_package(Catch2 3 REQUIRED)
  include(CTest)
  include(Catch)
  add_subdirectory(tests)
endif()

set(CMAKE_CXX_STANDARD 17)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -fbounds-check")
endif()

add_subdirectory(libspaecies)
add_subdirectory(rainshaft)
