find_package(Matplus++ 1.2.2 QUIET)
if (NOT Matplus++_FOUND)
  include(FetchContent)
  set(FETCHCONTENT_QUIET OFF)
  set(FETCHCONTENT_UPDATES_DISCONNECTED OFF)
  FetchContent_Declare(matplot++
    GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus.git
    GIT_TAG v1.2.2
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
  )
  set(MATPLOTPP_BUILD_SHARED_LIBS OFF)
  set(DMATPLOTPP_BUILD_TESTS OFF)
  FetchContent_MakeAvailable(matplot++)
endif()

add_executable(analyze_timescales analyze_timescales.cpp)

target_link_libraries(analyze_timescales
  rainshaft_common
  NetCDF::NetCDF
  SUNDIALS::arkode
  Boost::boost
  Boost::program_options
  spaecies
  Matplot++::matplot
)
target_include_directories(analyze_timescales PUBLIC ${CMAKE_SOURCE_DIR})
