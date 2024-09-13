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
  NetcdfWriter(const std::string& file_name);
  // Destructor closes the file.
  ~NetcdfWriter();

  // Write RainshaftGrid information to file.
  void write_grid(const RainshaftGrid& grid);

  // Write a series of States to file.
  void write_states(const std::vector<spaecies::State<const double>>& arrays);

  // Write a series of RainshaftDerivedVars to file.
  void write_derived_vars(const std::vector<RainshaftDerivedVars>& dvars);

  // Write total number of RHS evaluations to file.
  void write_num_rhs_evals(long int num_rhs_evals);

  // Write wallclock time for solution in ms to file.
  void write_walltime_ms(double walltime_ms);

protected:

  // NetCDF file id
  int ncid;

};

#endif // RAINSHAFT_NCIO_HPP
