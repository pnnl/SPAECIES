add_executable(rainshaft
  rainshaft_driver.cpp
  rainshaft_derived_vars.cpp
  rainshaft_integrator.cpp
  rainshaft_grid.cpp
  rainshaft_ncio.cpp
  rainshaft_state.cpp
  rainshaft_tendency.cpp
  sedimentation.cpp
)
target_link_libraries(rainshaft NetCDF::NetCDF spaecies)
target_include_directories(spaecies PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)
