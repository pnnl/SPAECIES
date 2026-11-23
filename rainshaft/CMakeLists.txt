add_executable(rainshaft
  evaporation.cpp
  explicit_integrator.cpp
  imex_integrator.cpp
  mri_integrator.cpp
  fixed_substep_integrator.cpp
  forcing_integrator.cpp
  limiting_integrator.cpp
  operator_splitting_integrator.cpp
  lookup_linear.cpp
  nudging.cpp
  nvector_serial_diagnostic.cpp
  parsing_utilities.cpp
  rainshaft_driver.cpp
  rainshaft_derived_vars.cpp
  rainshaft_grid.cpp
  rainshaft_initial.cpp
  rainshaft_ncio.cpp
  rainshaft_solution.cpp
  rainshaft_sum_process.cpp
  rainshaft_tend_descs.cpp
  saturation.cpp
  sed_cfl_integrator.cpp
  sedimentation.cpp
  self_collision.cpp
  sequential_split_integrator.cpp
  size_limiters.cpp
)
target_link_libraries(rainshaft
  NetCDF::NetCDF
  SUNDIALS::arkode
  SUNDIALS::sunlinsollapackdense
  Boost::boost
  Boost::program_options
  stdc++fs
  spaecies)
target_include_directories(rainshaft PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# CJV: can consider always creating the rainshaft_common object library (before
# the rainshaft executable) and having the rainshaft executable link against it
# to avoid recompiling the common source files when tools are enabled
if (BUILD_TOOLS)
  add_library(rainshaft_common OBJECT
    evaporation.cpp
    lookup_linear.cpp
    rainshaft_derived_vars.cpp
    rainshaft_grid.cpp
    rainshaft_tend_descs.cpp
    saturation.cpp
    sedimentation.cpp
    self_collision.cpp
  )
  target_include_directories(rainshaft_common PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_link_libraries(rainshaft_common
    SUNDIALS::arkode
    spaecies
  )
endif()
