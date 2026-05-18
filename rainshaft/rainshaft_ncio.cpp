#include "rainshaft_ncio.hpp"
#include <cstddef>
#include <vector>

NetcdfReader::NetcdfReader(const std::string& file_name) : is_open(true) {
  // SPS: Should throw error if this fails, probably also allow the file mode
  // argument to be adjusted.
  nc_open(file_name.c_str(), NC_NOWRITE, &fid);
}

NetcdfReader::~NetcdfReader() {
  // SPS: Can this fail? Probably safe to not check error code?
  if (is_open) {
    nc_close(fid);
  }
  is_open = false;
}

std::tuple<std::size_t, std::size_t> NetcdfReader::read_num_cases_and_max_levs() {
  int caseid, levid;
  nc_inq_dimid(fid, "nsamp", &caseid);
  nc_inq_dimid(fid, "nlev", &levid);
  std::size_t num_cases, max_levs;
  nc_inq_dimlen(fid, caseid, &num_cases);
  nc_inq_dimlen(fid, levid, &max_levs);
  return {num_cases, max_levs};
}

RainshaftGrid NetcdfReader::read_grid(std::size_t case_idx) {
  int nlevid, psid, pdelid;
  nc_inq_varid(fid, "nlev_samp", &nlevid);
  nc_inq_varid(fid, "surface_pressure", &psid);
  nc_inq_varid(fid, "pdel", &pdelid);
  int nlev;
  nc_get_var1_int(fid, nlevid, &case_idx, &nlev);
  std::size_t starts[2] = {case_idx, 0};
  std::size_t counts[2] = {1, (std::size_t) nlev};
  std::vector<double> pdel(nlev);
  nc_get_vara_double(fid, pdelid, starts, counts, pdel.data());
  // Dimension is nlev because we discard the top (cloud-base) level.
  std::vector<double> p_int(nlev);
  nc_get_var1_double(fid, psid, &case_idx, &p_int[nlev-1]);
  for (std::size_t i = nlev-1; i != 0; --i) {
    p_int[i-1] = p_int[i] - pdel[i];
  }
  return RainshaftGrid(p_int);
}

void NetcdfReader::read_boundary_conditions(std::size_t case_idx,
                                            RainshaftConstants &constants) {
  // Get state at cloud base.
  int rainfracid, pmidid, tid, qid, ncid, qcid, nrid, qrid;
  nc_inq_varid(fid, "rainfrac", &rainfracid);
  nc_inq_varid(fid, "pmid", &pmidid);
  nc_inq_varid(fid, "t", &tid);
  nc_inq_varid(fid, "q", &qid);
  nc_inq_varid(fid, "nc", &ncid);
  nc_inq_varid(fid, "qc", &qcid);
  nc_inq_varid(fid, "nr", &nrid);
  nc_inq_varid(fid, "qr", &qrid);
  double rainfrac, pmid, t, q, nc, qc, nr, qr;
  std::size_t loc[2] = {case_idx, 0};
  nc_get_var1_double(fid, rainfracid, loc, &rainfrac);
  nc_get_var1_double(fid, pmidid, loc, &pmid);
  nc_get_var1_double(fid, tid, loc, &t);
  nc_get_var1_double(fid, qid, loc, &q);
  nc_get_var1_double(fid, ncid, loc, &nc);
  nc_get_var1_double(fid, qcid, loc, &qc);
  nc_get_var1_double(fid, nrid, loc, &nr);
  nc_get_var1_double(fid, qrid, loc, &qr);
  constants.rho_top = rho_dry_from_ideal_gas_law(constants.rdry, constants.epsilon_h2o,
                                                 pmid, t, q);
  // Below assumes that the cloud fraction is 1.
  constants.nc_top = nc;
  constants.qc_top = qc;
  // Convert nr/qr at the top to in-cloud values.
  constants.nr_top = nr / rainfrac;
  constants.qr_top = qr / rainfrac;
}

void NetcdfReader::read_initial_conditions(std::size_t case_idx, State &initial_state) {
  int nlevid, rainfracid, tid, qid, ncid, qcid, nrid, qrid;
  nc_inq_varid(fid, "nlev_samp", &nlevid);
  nc_inq_varid(fid, "rainfrac", &rainfracid);
  nc_inq_varid(fid, "t", &tid);
  nc_inq_varid(fid, "q", &qid);
  nc_inq_varid(fid, "nc", &ncid);
  nc_inq_varid(fid, "qc", &qcid);
  nc_inq_varid(fid, "nr", &nrid);
  nc_inq_varid(fid, "qr", &qrid);
  int nlev;
  nc_get_var1_int(fid, nlevid, &case_idx, &nlev);
  std::vector<double> rainfrac(nlev-1);
  std::size_t starts[2] = {case_idx, 1};
  std::size_t counts[2] = {1, (std::size_t) nlev-1};
  nc_get_vara_double(fid, rainfracid, starts, counts, rainfrac.data());
  VarMut t = initial_state.get_variable("T");
  VarMut q = initial_state.get_variable("q");
  VarMut nc = initial_state.get_variable("nc");
  VarMut qc = initial_state.get_variable("qc");
  VarMut nr = initial_state.get_variable("nr");
  VarMut qr = initial_state.get_variable("qr");
  nc_get_vara_double(fid, tid, starts, counts, &t[0]);
  nc_get_vara_double(fid, qid, starts, counts, &q[0]);
  nc_get_vara_double(fid, ncid, starts, counts, &nc[0]);
  nc_get_vara_double(fid, qcid, starts, counts, &qc[0]);
  nc_get_vara_double(fid, nrid, starts, counts, &nr[0]);
  nc_get_vara_double(fid, qrid, starts, counts, &qr[0]);
  // Convert to in-cloud values.
  for (int i = 0; i != nlev - 1; ++i) {
    // Cloud fraction not implemented, but rain fraction is.
    nr[i] /= rainfrac[i];
    qr[i] /= rainfrac[i];
  }
}

NetcdfWriter::NetcdfWriter(const std::string& file_name, std::size_t num_cases, std::size_t max_levs) {
  // SPS: Should throw error if this fails, probably also allow the file mode
  // argument to be adjusted.
  nc_create(file_name.c_str(), NC_CLOBBER|NC_NETCDF4, &fid);
  int throwaway;
  nc_def_dim(fid, "case", num_cases, &throwaway);
  nc_def_dim(fid, "lev", max_levs, &throwaway);
  nc_def_dim(fid, "ilev", max_levs+1, &throwaway);
}

NetcdfWriter::~NetcdfWriter() {
  // SPS: Can this fail? Probably safe to not check error code?
  nc_close(fid);
}

void NetcdfWriter::write_grid(const RainshaftGrid& grid, std::size_t case_idx) {
  int caseid, levid, ilevid, nlevid, p_intid, p_midid;
  // SPS: Need to check errors from all of these calls.
  // SPS: Also should add metadata to variables being output.
  nc_inq_dimid(fid, "case", &caseid);
  nc_inq_dimid(fid, "lev", &levid);
  nc_inq_dimid(fid, "ilev", &ilevid);
  // Define variables.
  int status;
  status = nc_inq_varid(fid, "nlev", &nlevid);
  if (status == NC_ENOTVAR) {
    nc_def_var(fid, "nlev", NC_INT, 1, &caseid, &nlevid);
  }
  status = nc_inq_varid(fid, "p_int", &p_intid);
  if (status == NC_ENOTVAR) {
    int p_int_dimids[2] = {caseid, ilevid};
    nc_def_var(fid, "p_int", NC_DOUBLE, 2, p_int_dimids, &p_intid);
  }
  status = nc_inq_varid(fid, "p_mid", &p_midid);
  if (status == NC_ENOTVAR) {
    int p_mid_dimids[2] = {caseid, levid};
    nc_def_var(fid, "p_mid", NC_DOUBLE, 2, p_mid_dimids, &p_midid);
  }
  // Write variables.
  int nlev_int = grid.nlev;
  nc_put_var1_int(fid, nlevid, &case_idx, &nlev_int);
  std::size_t starts[2] = {case_idx, 0};
  std::size_t interface_counts[2] = {1, grid.nlev+1};
  std::size_t counts[2] = {1, grid.nlev};
  nc_put_vara_double(fid, p_intid, starts, interface_counts, grid.p_int.data());
  nc_put_vara_double(fid, p_midid, starts, counts, grid.p_mid.data());
}

void NetcdfWriter::write_states(const std::vector<StateConst>& arrays, std::size_t case_idx) {
  int caseid, levid, timeid;
  std::size_t nlev = arrays[0].get_variable("T").size();
  std::size_t ntimes = arrays.size();
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  nc_inq_dimid(fid, "case", &caseid);
  // SPS: Error if vertical coordinate not defined?
  // Find the vertical coordinate ids.
  nc_inq_dimid(fid, "lev", &levid);
  // SPS: Test this with and without a time dimension already present.
  int status = nc_inq_dimid(fid, "time", &timeid);
  if (status == NC_EBADDIM) {
    // Define time dimension.
    nc_def_dim(fid, "time", ntimes, &timeid);
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
    status = nc_inq_varid(fid, name.c_str(), &varid);
    if (status == NC_ENOTVAR) {
      nc_def_var(fid, name.c_str(), NC_DOUBLE, 3, var_dimids, &varid);
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
      nc_put_vara_double(fid, varids[j], starts, counts,
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
  nc_inq_dimid(fid, "case", &caseid);
  nc_inq_dimid(fid, "lev", &levid);
  nc_inq_dimid(fid, "ilev", &ilevid);
  nc_inq_dimid(fid, "time", &timeid);
  // Define derived variables.
  int status;
  int dimids[3] = {caseid, timeid, levid};
  int interface_dimids[3] = {caseid, timeid, ilevid};
  status = nc_inq_varid(fid, "z_int", &z_intid);
  if (status == NC_ENOTVAR) {
    nc_def_var(fid, "z_int", NC_DOUBLE, 3, interface_dimids, &z_intid);
  }
  status = nc_inq_varid(fid, "dz", &dzid);
  if (status == NC_ENOTVAR) {
    nc_def_var(fid, "dz", NC_DOUBLE, 3, dimids, &dzid);
  }
  status = nc_inq_varid(fid, "rho_dry", &rho_dryid);
  if (status == NC_ENOTVAR) {
    nc_def_var(fid, "rho_dry", NC_DOUBLE, 3, dimids, &rho_dryid);
  }
  status = nc_inq_varid(fid, "lambdar", &lambdarid);
  if (status == NC_ENOTVAR) {
    nc_def_var(fid, "lambdar", NC_DOUBLE, 3, dimids, &lambdarid);
  }

  // Write variables.
  for (std::size_t i = 0; i != ntimes; ++i) {
    std::size_t starts[3] = {case_idx, i, 0};
    std::size_t counts[3] = {1, 1, nlev};
    std::size_t interface_counts[3] = {1, 1, nlev+1};
    nc_put_vara_double(fid, z_intid, starts, interface_counts, dvars[i].z_int.data());
    nc_put_vara_double(fid, dzid, starts, counts, dvars[i].dz.data());
    nc_put_vara_double(fid, rho_dryid, starts, counts, dvars[i].rho_dry.data());
    nc_put_vara_double(fid, lambdarid, starts, counts, dvars[i].lambdar.data());
  }
}

void NetcdfWriter::write_num_rhs_evals(std::int64_t num_rhs_evals, std::size_t case_idx) {
  int caseid, evalsid;
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  nc_inq_dimid(fid, "case", &caseid);
  // Define derived variables.
  int status = nc_inq_varid(fid, "num_rhs_evals", &evalsid);
  if (status == NC_ENOTVAR) {
    nc_def_var(fid, "num_rhs_evals", NC_INT64, 1, &caseid, &evalsid);
  }
  // Write variable.
  // Use the generic output rather than rely on "long" or "longlong" being
  // a particular length.
  nc_put_var1(fid, evalsid, &case_idx, &num_rhs_evals);
}

void NetcdfWriter::write_walltime_ms(double walltime_ms, std::size_t case_idx) {
  int caseid, walltime_msid;
  // SPS: Need to check errors from all these as well.
  // SPS: Add variable metadata to all these too.
  nc_inq_dimid(fid, "case", &caseid);
  // Define derived variables.
  int status = nc_inq_varid(fid, "walltime_ms", &walltime_msid);
  if (status == NC_ENOTVAR) {
    nc_def_var(fid, "walltime_ms", NC_DOUBLE, 1, &caseid, &walltime_msid);
  }
  // Write variable.
  nc_put_var1_double(fid, walltime_msid, &case_idx, &walltime_ms);
}

void NetcdfWriter::write_metadata(int order, double dt, double dt_partition_1, double dt_partition_2, double rel_tol,
                                  bool postprocess, bool use_lookup, std::string method_type, int steps_per_output,
                                  std::string initial_condition_file, int num_cases, int icase_in, double final_time,
                                  bool do_nudging) {
  int orderid, dtid, dt_partition_1id, dt_partition_2id, rel_tolid;
  int postprocessid, use_lookupid, method_typeid, steps_per_outputid;
  int initial_condition_fileid, num_casesid, icase_inid, final_timeid, do_nudgingid;

  // method order
  nc_def_var(fid, "method_order", NC_INT, 0, NULL, &orderid);
  nc_put_var_int(fid, orderid, &order);
  
  // coupling time step
  nc_def_var(fid, "dt", NC_DOUBLE, 0, NULL, &dtid);
  nc_put_var_double(fid, dtid, &dt);

  // partition time steps
  nc_def_var(fid, "dt_partition_1", NC_DOUBLE, 0, NULL, &dt_partition_1id);
  nc_put_var_double(fid, dt_partition_1id, &dt_partition_1);
  nc_def_var(fid, "dt_partition_2", NC_DOUBLE, 0, NULL, &dt_partition_2id);
  nc_put_var_double(fid, dt_partition_2id, &dt_partition_2);

  // rel tol
  nc_def_var(fid, "rel_tol", NC_DOUBLE, 0, NULL, &rel_tolid);
  nc_put_var_double(fid, rel_tolid, &rel_tol);

  // limiter and lookup tables
  nc_def_var(fid, "postprocess", NC_UBYTE, 0, NULL, &postprocessid);
  nc_put_var(fid, postprocessid, &postprocess);
  nc_def_var(fid, "use_lookup", NC_UBYTE, 0, NULL, &use_lookupid);
  nc_put_var(fid, use_lookupid, &use_lookup);

  // type of integrator
  nc_def_var(fid, "method_type", NC_STRING, 0, NULL, &method_typeid);
  { // Logic to deal with the NC_STRING interface expecting char**.
    char* temp = new char[method_type.length()+1];
    std::strcpy(temp, method_type.c_str());
    nc_put_var(fid, method_typeid, &temp);
    delete[] temp;
  }

  // frequency of solution states
  nc_def_var(fid, "steps_per_output", NC_INT, 0, NULL, &steps_per_outputid);
  nc_put_var_int(fid, steps_per_outputid, &steps_per_output);

  // initial condition file
  nc_def_var(fid, "initial_condition_file", NC_STRING, 0, NULL, &initial_condition_fileid);
  { // Logic to deal with the NC_STRING interface expecting char**.
    char* temp = new char[initial_condition_file.length()+1];
    std::strcpy(temp, initial_condition_file.c_str());
    nc_put_var(fid, initial_condition_fileid, &temp);
    delete[] temp;
  }

  // number of cases simulated in this file
  nc_def_var(fid, "num_cases", NC_INT, 0, NULL, &num_casesid);
  nc_put_var_int(fid, num_casesid, &num_cases);

  // case index (-1 for all columns)
  nc_def_var(fid, "icase_in", NC_INT, 0, NULL, &icase_inid);
  nc_put_var_int(fid, icase_inid, &icase_in);

  // final time
  nc_def_var(fid, "final_time", NC_DOUBLE, 0, NULL, &final_timeid);
  nc_put_var_double(fid, final_timeid, &final_time);

  // nudging flag
  nc_def_var(fid, "do_nudging", NC_UBYTE, 0, NULL, &do_nudgingid);
  nc_put_var(fid, do_nudgingid, &do_nudging);
}
