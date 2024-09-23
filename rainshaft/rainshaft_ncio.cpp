#include "rainshaft_ncio.hpp"
#include <cstddef>

NetcdfWriter::NetcdfWriter(const std::string& file_name) {
  // SPS: Should throw error if this fails, probably also allow the file mode
  // argument to be adjusted.
  nc_create(file_name.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);
}

NetcdfWriter::~NetcdfWriter() {
  // SPS: Can this fail? Probably safe to not check error code?
  nc_close(ncid);
}

void NetcdfWriter::write_grid(const RainshaftGrid& grid) {
  int levid, ilevid, p_intid, p_midid;
  // SPS: Need to check errors from all of these calls.
  // SPS: Also should add metadata to variables being output.
  // Define dimensions.
  nc_def_dim(ncid, "lev", grid.nlev, &levid);
  nc_def_dim(ncid, "ilev", grid.nlev + 1, &ilevid);
  // Define variables.
  nc_def_var(ncid, "p_int", NC_DOUBLE, 1, &ilevid, &p_intid);
  nc_def_var(ncid, "p_mid", NC_DOUBLE, 1, &levid, &p_midid);

  // Write variables.
  nc_put_var(ncid, p_intid, grid.p_int.data());
  nc_put_var(ncid, p_midid, grid.p_mid.data());
}

void NetcdfWriter::write_states(const std::vector<StateConst>& arrays) {
  int levid, timeid;
  std::size_t nlev = arrays[0].get_variable("T").size();
  std::size_t ntimes = arrays.size();
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  // SPS: Error if vertical coordinate not defined?
  // Find the vertical coordinate ids.
  nc_inq_dimid(ncid, "lev", &levid);
  // SPS: Test this with and without a time dimension already present.
  int status = nc_inq_dimid(ncid, "time", &timeid);
  if (status == NC_EBADDIM) {
    // Define time dimension.
    nc_def_dim(ncid, "time", ntimes, &timeid);
  }
  // Define variables.
  // SPS: Should be more flexible about dimensions here.
  int var_dimids[2] = {timeid, levid};
  std::vector<std::string> varnames;
  std::vector<int> varids;
  // SPS: Need to check that all arrays have the same variables, or come up with
  // a new type for time series of arrays.
  for (spaecies::VarDescPtr var_desc : arrays[0].var_descs()) {
    int varid;
    std::string name = var_desc->name;
    nc_def_var(ncid, name.c_str(), NC_DOUBLE, 2, var_dimids, &varid);
    varnames.push_back(name);
    varids.push_back(varid);
  }

  // Write variables.
  for (std::size_t i = 0; i != ntimes; ++i) {
    std::size_t starts[2] = {i, 0};
    std::size_t counts[2] = {1, nlev};
    for (std::size_t j = 0; j != varids.size(); ++j) {
      // SPS: Would be more efficient/simple if we could request data using an
      // integer id from the VariableArrayView, rather than using string lookups.
      // SPS: Note also the implicit assumption that the variable data is
      // contiguous.
      auto var = arrays[i].get_variable(varnames[j]);
      nc_put_vara_double(ncid, varids[j], starts, counts,
                         &var[0]);
    }
  }
}

void NetcdfWriter::write_derived_vars(const std::vector<RainshaftDerivedVars>& dvars) {
  int levid, ilevid, timeid, z_intid, dzid, rho_dryid, lambdarid;
  std::size_t nlev = dvars[0].dz.size();
  std::size_t ntimes = dvars.size();
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  // SPS: Error if time and vertical coordinates not defined?
  // Find the vertical and time coordinate ids.
  nc_inq_dimid(ncid, "lev", &levid);
  nc_inq_dimid(ncid, "ilev", &ilevid);
  nc_inq_dimid(ncid, "time", &timeid);
  // Define derived variables.
  int dimids[2] = {timeid, levid};
  int interface_dimids[2] = {timeid, ilevid};
  nc_def_var(ncid, "z_int", NC_DOUBLE, 2, interface_dimids, &z_intid);
  nc_def_var(ncid, "dz", NC_DOUBLE, 2, dimids, &dzid);
  nc_def_var(ncid, "rho_dry", NC_DOUBLE, 2, dimids, &rho_dryid);
  nc_def_var(ncid, "lambdar", NC_DOUBLE, 2, dimids, &lambdarid);

  // Write variables.
  for (std::size_t i = 0; i != ntimes; ++i) {
    std::size_t starts[2] = {i, 0};
    std::size_t counts[2] = {1, nlev};
    std::size_t interface_counts[2] = {1, nlev+1};
    nc_put_vara_double(ncid, z_intid, starts, interface_counts, dvars[i].z_int.data());
    nc_put_vara_double(ncid, dzid, starts, counts, dvars[i].dz.data());
    nc_put_vara_double(ncid, rho_dryid, starts, counts, dvars[i].rho_dry.data());
    nc_put_vara_double(ncid, lambdarid, starts, counts, dvars[i].lambdar.data());
  }
}

void NetcdfWriter::write_num_rhs_evals(long int num_rhs_evals) {
  int evalsid;
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  // Define derived variables.
  nc_def_var(ncid, "num_rhs_evals", NC_INT, 0, NULL, &evalsid);

  // Write variable.
  nc_put_var(ncid, evalsid, &num_rhs_evals);
}

void NetcdfWriter::write_walltime_ms(double walltime_ms) {
  int walltime_msid;
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  // Define derived variables.
  nc_def_var(ncid, "walltime_ms", NC_DOUBLE, 0, NULL, &walltime_msid);

  // Write variable.
  nc_put_var(ncid, walltime_msid, &walltime_ms);
}
