cmake_minimum_required(VERSION 3.22)
project(SPAECIES LANGUAGES CXX VERSION 0.0.1)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(NetCDF REQUIRED)
find_package(SUNDIALS REQUIRED)
find_package(Boost REQUIRED)
find_package(Catch2 2...<3 REQUIRED)

include(CTest)
include(Catch)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_DEBUG "-Wall -fbounds-check")
endif()

add_subdirectory(libspaecies)
add_subdirectory(rainshaft)
add_subdirectory(tests)
