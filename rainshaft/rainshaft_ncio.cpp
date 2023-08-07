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
  int levid, ilevid, p_intid, p_midid, z_intid, dzid;
  // SPS: Need to check errors from all of these calls.
  // SPS: Also should add metadata to variables being output.
  // Define dimensions.
  nc_def_dim(ncid, "lev", grid.nlev, &levid);
  nc_def_dim(ncid, "ilev", grid.nlev + 1, &ilevid);
  // Define variables.
  nc_def_var(ncid, "p_int", NC_DOUBLE, 1, &ilevid, &p_intid);
  nc_def_var(ncid, "p_mid", NC_DOUBLE, 1, &levid, &p_midid);
  nc_def_var(ncid, "z_int", NC_DOUBLE, 1, &ilevid, &z_intid);
  nc_def_var(ncid, "dz", NC_DOUBLE, 1, &levid, &dzid);

  // Write variables.
  nc_put_var(ncid, p_intid, grid.p_int.data());
  nc_put_var(ncid, p_midid, grid.p_mid.data());
  nc_put_var(ncid, z_intid, grid.z_int.data());
  nc_put_var(ncid, dzid, grid.dz.data());
}

void NetcdfWriter::write_states(const std::vector<RainshaftState>& states) {
  int levid, timeid, tid, qid, nrid, qrid;
  std::size_t nlevs = states[0].t.size();
  std::size_t ntimes = states.size();
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  // SPS: Error if vertical coordinate not defined?
  // Find the vertical coordinate id.
  nc_inq_dimid(ncid, "lev", &levid);
  // Define time dimension.
  nc_def_dim(ncid, "time", ntimes, &timeid);
  // Define state variables.
  int state_dimids[2] = {timeid, levid};
  nc_def_var(ncid, "T", NC_DOUBLE, 2, state_dimids, &tid);
  nc_def_var(ncid, "q", NC_DOUBLE, 2, state_dimids, &qid);
  nc_def_var(ncid, "nr", NC_DOUBLE, 2, state_dimids, &nrid);
  nc_def_var(ncid, "qr", NC_DOUBLE, 2, state_dimids, &qrid);

  // Write variables.
  for (std::size_t i = 0; i != ntimes; ++i) {
    std::size_t starts[2] = {i, 0};
    std::size_t counts[2] = {1, nlevs};
    nc_put_vara_double(ncid, tid, starts, counts, states[i].t.data());
    nc_put_vara_double(ncid, qid, starts, counts, states[i].q.data());
    nc_put_vara_double(ncid, nrid, starts, counts, states[i].nr.data());
    nc_put_vara_double(ncid, qrid, starts, counts, states[i].qr.data());
  }
}

void NetcdfWriter::write_derived_vars(const std::vector<RainshaftDerivedVars>& dvars) {
  int levid, timeid, rhoid, lambdarid;
  std::size_t nlevs = dvars[0].rho.size();
  std::size_t ntimes = dvars.size();
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  // SPS: Error if time and vertical coordinates not defined?
  // Find the vertical and time coordinate ids.
  nc_inq_dimid(ncid, "lev", &levid);
  nc_inq_dimid(ncid, "time", &timeid);
  // Define derived variables.
  int dimids[2] = {timeid, levid};
  nc_def_var(ncid, "rho", NC_DOUBLE, 2, dimids, &rhoid);
  nc_def_var(ncid, "lambdar", NC_DOUBLE, 2, dimids, &lambdarid);

  // Write variables.
  for (std::size_t i = 0; i != ntimes; ++i) {
    std::size_t starts[2] = {i, 0};
    std::size_t counts[2] = {1, nlevs};
    nc_put_vara_double(ncid, rhoid, starts, counts, dvars[i].rho.data());
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
