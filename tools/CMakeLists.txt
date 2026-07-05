add_executable(analyze_timescales analyze_timescales.cpp)

target_link_libraries(analyze_timescales
  rainshaft_common
  NetCDF::NetCDF
  SUNDIALS::arkode
  Boost::boost
  Boost::program_options
  spaecies
)
target_include_directories(analyze_timescales PUBLIC ${CMAKE_SOURCE_DIR})
