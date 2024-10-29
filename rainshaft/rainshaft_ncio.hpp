#ifndef RAINSHAFT_NCIO_HPP
#define RAINSHAFT_NCIO_HPP
#include <string>
#include <vector>

extern "C"
{
#include "netcdf.h"
}

#include "spaecies.hpp"

#include "rainshaft_grid.hpp"
#include "rainshaft_derived_vars.hpp"

class NetcdfWriter {
  
public:

  // Constructed from file name.
  NetcdfWriter(const std::string& file_name, std::size_t num_cases, std::size_t max_levs);
  // Destructor closes the file.
  ~NetcdfWriter();

  // Write RainshaftGrid information to file.
  void write_grid(const RainshaftGrid& grid, std::size_t case_idx);

  // Write a series of States to file.
  void write_states(const std::vector<StateConst>& arrays, std::size_t case_idx);

  // Write a series of RainshaftDerivedVars to file.
  void write_derived_vars(const std::vector<RainshaftDerivedVars>& dvars, std::size_t case_idx);

  // Write total number of RHS evaluations to file.
  void write_num_rhs_evals(std::int64_t num_rhs_evals, std::size_t case_idx);

  // Write wallclock time for solution in ms to file.
  void write_walltime_ms(double walltime_ms, std::size_t case_idx);

protected:

  // NetCDF file id
  int ncid;

};

#endif // RAINSHAFT_NCIO_HPP
