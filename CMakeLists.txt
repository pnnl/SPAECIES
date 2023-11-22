cmake_minimum_required(VERSION 3.22)

if( "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" )
  message( FATAL_ERROR
    " FATAL: In-source builds are not allowed."
    " Remember to delete CMakeCache.txt and any other generated files"
    " (e.g. the CMakeFiles/ directory), then build in a separate source tree.")
endif()

project(SPAECIES LANGUAGES CXX VERSION 0.0.1)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(NetCDF REQUIRED)
find_package(SUNDIALS 7.1.0 REQUIRED)
find_package(Boost REQUIRED)

if (BUILD_TESTING)
  find_package(Catch2 3 REQUIRED)
  include(CTest)
  include(Catch)
  add_subdirectory(tests)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS_DEBUG "-Wall -fbounds-check -g")
  set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif()

add_subdirectory(libspaecies)
add_subdirectory(rainshaft)
