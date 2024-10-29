#include "rainshaft_ncio.hpp"
#include <cstddef>

NetcdfWriter::NetcdfWriter(const std::string& file_name, std::size_t num_cases, std::size_t max_levs) {
  // SPS: Should throw error if this fails, probably also allow the file mode
  // argument to be adjusted.
  nc_create(file_name.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);
  int throwaway;
  nc_def_dim(ncid, "case", num_cases, &throwaway);
  nc_def_dim(ncid, "lev", max_levs, &throwaway);
  nc_def_dim(ncid, "ilev", max_levs+1, &throwaway);
}

NetcdfWriter::~NetcdfWriter() {
  // SPS: Can this fail? Probably safe to not check error code?
  nc_close(ncid);
}

void NetcdfWriter::write_grid(const RainshaftGrid& grid, std::size_t case_idx) {
  int caseid, levid, ilevid, nlevid, p_intid, p_midid;
  // SPS: Need to check errors from all of these calls.
  // SPS: Also should add metadata to variables being output.
  nc_inq_dimid(ncid, "case", &caseid);
  nc_inq_dimid(ncid, "lev", &levid);
  nc_inq_dimid(ncid, "ilev", &ilevid);
  // Define variables.
  int status;
  status = nc_inq_varid(ncid, "nlev", &nlevid);
  if (status == NC_ENOTVAR) {
    nc_def_var(ncid, "nlev", NC_INT, 1, &caseid, &nlevid);
  }
  status = nc_inq_varid(ncid, "p_int", &p_intid);
  if (status == NC_ENOTVAR) {
    int p_int_dimids[2] = {caseid, ilevid};
    nc_def_var(ncid, "p_int", NC_DOUBLE, 2, p_int_dimids, &p_intid);
  }
  int p_int_dimids[2] = {caseid, ilevid};
  status = nc_inq_varid(ncid, "p_mid", &p_midid);
  if (status == NC_ENOTVAR) {
    int p_mid_dimids[2] = {caseid, levid};
    nc_def_var(ncid, "p_mid", NC_DOUBLE, 2, p_mid_dimids, &p_midid);
  }
  // Write variables.
  int nlev_int = grid.nlev;
  nc_put_var1_int(ncid, nlevid, &case_idx, &nlev_int);
  std::size_t starts[2] = {case_idx, 0};
  std::size_t interface_counts[2] = {1, grid.nlev+1};
  std::size_t counts[2] = {1, grid.nlev};
  nc_put_vara_double(ncid, p_intid, starts, interface_counts, grid.p_int.data());
  nc_put_vara_double(ncid, p_midid, starts, counts, grid.p_mid.data());
}

void NetcdfWriter::write_states(const std::vector<StateConst>& arrays, std::size_t case_idx) {
  int caseid, levid, timeid;
  std::size_t nlev = arrays[0].get_variable("T").size();
  std::size_t ntimes = arrays.size();
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  nc_inq_dimid(ncid, "case", &caseid);
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
  int var_dimids[3] = {caseid, timeid, levid};
  std::vector<std::string> varnames;
  std::vector<int> varids;
  // SPS: Need to check that all arrays have the same variables, or come up with
  // a new type for time series of arrays.
  for (spaecies::VarDescPtr var_desc : arrays[0].var_descs()) {
    int varid;
    std::string name = var_desc->name;
    status = nc_inq_varid(ncid, name.c_str(), &varid);
    if (status == NC_ENOTVAR) {
      nc_def_var(ncid, name.c_str(), NC_DOUBLE, 3, var_dimids, &varid);
    }
    varnames.push_back(name);
    varids.push_back(varid);
  }

  // Write variables.
  for (std::size_t i = 0; i != ntimes; ++i) {
    std::size_t starts[3] = {case_idx, i, 0};
    std::size_t counts[3] = {1, 1, nlev};
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

void NetcdfWriter::write_derived_vars(const std::vector<RainshaftDerivedVars>& dvars, std::size_t case_idx) {
  int caseid, levid, ilevid, timeid, z_intid, dzid, rho_dryid, lambdarid;
  std::size_t nlev = dvars[0].dz.size();
  std::size_t ntimes = dvars.size();
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  // SPS: Error if time and vertical coordinates not defined?
  // Find the vertical and time coordinate ids.
  nc_inq_dimid(ncid, "case", &caseid);
  nc_inq_dimid(ncid, "lev", &levid);
  nc_inq_dimid(ncid, "ilev", &ilevid);
  nc_inq_dimid(ncid, "time", &timeid);
  // Define derived variables.
  int status;
  int dimids[3] = {caseid, timeid, levid};
  int interface_dimids[3] = {caseid, timeid, ilevid};
  status = nc_inq_varid(ncid, "z_int", &z_intid);
  if (status == NC_ENOTVAR) {
    nc_def_var(ncid, "z_int", NC_DOUBLE, 3, interface_dimids, &z_intid);
  }
  status = nc_inq_varid(ncid, "dz", &dzid);
  if (status == NC_ENOTVAR) {
    nc_def_var(ncid, "dz", NC_DOUBLE, 3, dimids, &dzid);
  }
  status = nc_inq_varid(ncid, "rho_dry", &rho_dryid);
  if (status == NC_ENOTVAR) {
    nc_def_var(ncid, "rho_dry", NC_DOUBLE, 3, dimids, &rho_dryid);
  }
  status = nc_inq_varid(ncid, "lambdar", &lambdarid);
  if (status == NC_ENOTVAR) {
    nc_def_var(ncid, "lambdar", NC_DOUBLE, 3, dimids, &lambdarid);
  }

  // Write variables.
  for (std::size_t i = 0; i != ntimes; ++i) {
    std::size_t starts[3] = {case_idx, i, 0};
    std::size_t counts[3] = {1, 1, nlev};
    std::size_t interface_counts[3] = {1, 1, nlev+1};
    nc_put_vara_double(ncid, z_intid, starts, interface_counts, dvars[i].z_int.data());
    nc_put_vara_double(ncid, dzid, starts, counts, dvars[i].dz.data());
    nc_put_vara_double(ncid, rho_dryid, starts, counts, dvars[i].rho_dry.data());
    nc_put_vara_double(ncid, lambdarid, starts, counts, dvars[i].lambdar.data());
  }
}

void NetcdfWriter::write_num_rhs_evals(std::int64_t num_rhs_evals, std::size_t case_idx) {
  int caseid, evalsid;
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  nc_inq_dimid(ncid, "case", &caseid);
  // Define derived variables.
  int status = nc_inq_varid(ncid, "num_rhs_evals", &evalsid);
  if (status == NC_ENOTVAR) {
    nc_def_var(ncid, "num_rhs_evals", NC_INT64, 1, &caseid, &evalsid);
  }
  // Write variable.
  nc_put_var1_longlong(ncid, evalsid, &case_idx, &num_rhs_evals);
}

void NetcdfWriter::write_walltime_ms(double walltime_ms, std::size_t case_idx) {
  int caseid, walltime_msid;
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  nc_inq_dimid(ncid, "case", &caseid);
  // Define derived variables.
  int status = nc_inq_varid(ncid, "walltime_ms", &walltime_msid);
  if (status == NC_ENOTVAR) {
    nc_def_var(ncid, "walltime_ms", NC_DOUBLE, 1, &caseid, &walltime_msid);
  }
  // Write variable.
  nc_put_var1_double(ncid, walltime_msid, &case_idx, &walltime_ms);
}
