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

#include <tuple>

class NetcdfReader {

public:

  // Constructed from file name.
  NetcdfReader(const std::string& file_name);
  // Destructor closes the file.
  ~NetcdfReader();
  // Other constructors as per the rule of 5.
  NetcdfReader(const NetcdfReader&) = delete;
  NetcdfReader& operator=(const NetcdfReader &) = delete;
  NetcdfReader(NetcdfReader &&other) noexcept : is_open(true), ncid(other.ncid) {
    other.is_open = false;
    other.ncid = -100;
  }
  NetcdfReader& operator=(NetcdfReader &&other) {
    is_open = true;
    ncid = other.ncid;
    other.is_open = false;
    other.ncid = -100;
    return *this;
  }

  std::tuple<std::size_t, std::size_t> read_num_cases_and_max_levs();

  RainshaftGrid read_grid(std::size_t case_idx);

  void read_boundary_conditions(std::size_t case_idx, RainshaftConstants &constants);

  void read_initial_conditions(std::size_t case_idx, State &initial_state);

protected:

  // Is the file open?
  bool is_open;

  // NetCDF file id
  int ncid;

};

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

  // Write all user inputs to file
  void write_metadata(int order, double dt, double dt_partition_1, double dt_partition_2, double rel_tol,
                      bool postprocess, bool use_lookup, std::string method_type, int steps_per_output,
                      std::string initial_condition_file, int num_cases, int icase_in, double final_time,
                      bool do_nudging);
protected:

  // NetCDF file id
  int ncid;

};

#endif // RAINSHAFT_NCIO_HPP
