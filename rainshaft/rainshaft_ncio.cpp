#include "rainshaft_ncio.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check_netcdf(const int status, const std::string& operation) {
  if (status != NC_NOERR) {
    throw std::runtime_error("NetCDF " + operation + " failed: " + nc_strerror(status));
  }
}

int get_dimension_id(const int ncid, const char* name) {
  int dimid;
  check_netcdf(nc_inq_dimid(ncid, name, &dimid),
               "query dimension '" + std::string(name) + "'");
  return dimid;
}

int get_variable_id(const int ncid, const char* name) {
  int varid;
  check_netcdf(nc_inq_varid(ncid, name, &varid),
               "query variable '" + std::string(name) + "'");
  return varid;
}

int get_or_define_variable(const int ncid, const char* name, const nc_type type,
                           const int ndims, const int* dimids) {
  int varid;
  const int status = nc_inq_varid(ncid, name, &varid);
  if (status == NC_ENOTVAR) {
    check_netcdf(nc_def_var(ncid, name, type, ndims, dimids, &varid),
                 "define variable '" + std::string(name) + "'");
  } else {
    check_netcdf(status, "query variable '" + std::string(name) + "'");
  }
  return varid;
}

int define_scalar_variable(const int ncid, const char* name, const nc_type type) {
  int varid;
  check_netcdf(nc_def_var(ncid, name, type, 0, nullptr, &varid),
               "define scalar variable '" + std::string(name) + "'");
  return varid;
}

void require_nonempty(const std::size_t size, const char* name) {
  if (size == 0) {
    throw std::invalid_argument(std::string(name) + " must not be empty");
  }
}

} // namespace

NetcdfReader::NetcdfReader(const std::string& file_name) {
  check_netcdf(nc_open(file_name.c_str(), NC_NOWRITE, &ncid),
               "open input file '" + file_name + "'");
}

NetcdfReader::~NetcdfReader() {
  close();
}

void NetcdfReader::close() noexcept {
  if (ncid != -1) {
    nc_close(ncid);
    ncid = -1;
  }
}

std::tuple<std::size_t, std::size_t> NetcdfReader::read_num_cases_and_max_levs() {
  const int caseid = get_dimension_id(ncid, "nsamp");
  const int levid = get_dimension_id(ncid, "nlev");
  std::size_t num_cases, max_levs;
  check_netcdf(nc_inq_dimlen(ncid, caseid, &num_cases), "query length of dimension 'nsamp'");
  check_netcdf(nc_inq_dimlen(ncid, levid, &max_levs), "query length of dimension 'nlev'");
  return {num_cases, max_levs};
}

RainshaftGrid NetcdfReader::read_grid(const std::size_t case_idx) {
  const int nlevid = get_variable_id(ncid, "nlev_samp");
  const int psid = get_variable_id(ncid, "surface_pressure");
  const int pdelid = get_variable_id(ncid, "pdel");
  int nlev;
  check_netcdf(nc_get_var1_int(ncid, nlevid, &case_idx, &nlev), "read variable 'nlev_samp'");
  if (nlev <= 0) {
    throw std::runtime_error("NetCDF variable 'nlev_samp' must be positive");
  }

  const std::size_t starts[2] = {case_idx, 0};
  const std::size_t counts[2] = {1, static_cast<std::size_t>(nlev)};
  std::vector<double> pdel(nlev);
  check_netcdf(nc_get_vara_double(ncid, pdelid, starts, counts, pdel.data()),
               "read variable 'pdel'");

  // The cloud-base level is discarded, leaving nlev model interfaces.
  std::vector<double> p_int(nlev);
  check_netcdf(nc_get_var1_double(ncid, psid, &case_idx, &p_int.back()),
               "read variable 'surface_pressure'");
  for (std::size_t i = p_int.size() - 1; i != 0; --i) {
    p_int[i - 1] = p_int[i] - pdel[i];
  }
  return RainshaftGrid(p_int);
}

void NetcdfReader::read_boundary_conditions(const std::size_t case_idx,
                                            RainshaftConstants& constants) {
  const int rainfracid = get_variable_id(ncid, "rainfrac");
  const int pmidid = get_variable_id(ncid, "pmid");
  const int tid = get_variable_id(ncid, "t");
  const int qid = get_variable_id(ncid, "q");
  const int nrid = get_variable_id(ncid, "nr");
  const int qrid = get_variable_id(ncid, "qr");
  double rainfrac, pmid, t, q, nr, qr;
  const std::size_t location[2] = {case_idx, 0};
  check_netcdf(nc_get_var1_double(ncid, rainfracid, location, &rainfrac), "read variable 'rainfrac'");
  check_netcdf(nc_get_var1_double(ncid, pmidid, location, &pmid), "read variable 'pmid'");
  check_netcdf(nc_get_var1_double(ncid, tid, location, &t), "read variable 't'");
  check_netcdf(nc_get_var1_double(ncid, qid, location, &q), "read variable 'q'");
  check_netcdf(nc_get_var1_double(ncid, nrid, location, &nr), "read variable 'nr'");
  check_netcdf(nc_get_var1_double(ncid, qrid, location, &qr), "read variable 'qr'");

  constants.rho_top = rho_dry_from_ideal_gas_law(constants.rdry, constants.epsilon_h2o,
                                                  pmid, t, q);
  constants.nr_top = nr / rainfrac;
  constants.qr_top = qr / rainfrac;
}

void NetcdfReader::read_initial_conditions(const std::size_t case_idx, State& initial_state) {
  const int nlevid = get_variable_id(ncid, "nlev_samp");
  const int rainfracid = get_variable_id(ncid, "rainfrac");
  const int tid = get_variable_id(ncid, "t");
  const int qid = get_variable_id(ncid, "q");
  const int nrid = get_variable_id(ncid, "nr");
  const int qrid = get_variable_id(ncid, "qr");
  int nlev;
  check_netcdf(nc_get_var1_int(ncid, nlevid, &case_idx, &nlev), "read variable 'nlev_samp'");
  if (nlev < 2) {
    throw std::runtime_error("NetCDF variable 'nlev_samp' must be at least 2");
  }

  const std::size_t starts[2] = {case_idx, 1};
  const std::size_t counts[2] = {1, static_cast<std::size_t>(nlev - 1)};
  std::vector<double> rainfrac(nlev - 1);
  check_netcdf(nc_get_vara_double(ncid, rainfracid, starts, counts, rainfrac.data()),
               "read variable 'rainfrac'");

  VarMut t = initial_state.get_variable("T").value();
  VarMut q = initial_state.get_variable("q").value();
  VarMut nr = initial_state.get_variable("nr").value();
  VarMut qr = initial_state.get_variable("qr").value();
  check_netcdf(nc_get_vara_double(ncid, tid, starts, counts, &t[0]), "read variable 't'");
  check_netcdf(nc_get_vara_double(ncid, qid, starts, counts, &q[0]), "read variable 'q'");
  check_netcdf(nc_get_vara_double(ncid, nrid, starts, counts, &nr[0]), "read variable 'nr'");
  check_netcdf(nc_get_vara_double(ncid, qrid, starts, counts, &qr[0]), "read variable 'qr'");

  for (std::size_t i = 0; i != rainfrac.size(); ++i) {
    nr[i] /= rainfrac[i];
    qr[i] /= rainfrac[i];
  }
}

NetcdfWriter::NetcdfWriter(const std::string& file_name, const std::size_t num_cases,
                           const std::size_t max_levs) {
  check_netcdf(nc_create(file_name.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid),
               "create output file '" + file_name + "'");
  int dimid;
  check_netcdf(nc_def_dim(ncid, "case", num_cases, &dimid), "define dimension 'case'");
  check_netcdf(nc_def_dim(ncid, "lev", max_levs, &dimid), "define dimension 'lev'");
  check_netcdf(nc_def_dim(ncid, "ilev", max_levs + 1, &dimid), "define dimension 'ilev'");
}

NetcdfWriter::~NetcdfWriter() {
  close();
}

void NetcdfWriter::close() noexcept {
  if (ncid != -1) {
    nc_close(ncid);
    ncid = -1;
  }
}

void NetcdfWriter::write_grid(const RainshaftGrid& grid, const std::size_t case_idx) {
  const int caseid = get_dimension_id(ncid, "case");
  const int levid = get_dimension_id(ncid, "lev");
  const int ilevid = get_dimension_id(ncid, "ilev");
  const int nlevid = get_or_define_variable(ncid, "nlev", NC_INT, 1, &caseid);
  const int p_int_dimids[2] = {caseid, ilevid};
  const int p_intid = get_or_define_variable(ncid, "p_int", NC_DOUBLE, 2, p_int_dimids);
  const int p_mid_dimids[2] = {caseid, levid};
  const int p_midid = get_or_define_variable(ncid, "p_mid", NC_DOUBLE, 2, p_mid_dimids);

  const int nlev = grid.nlev;
  check_netcdf(nc_put_var1_int(ncid, nlevid, &case_idx, &nlev), "write variable 'nlev'");
  const std::size_t starts[2] = {case_idx, 0};
  const std::size_t interface_counts[2] = {1, grid.nlev + 1};
  const std::size_t counts[2] = {1, grid.nlev};
  check_netcdf(nc_put_vara_double(ncid, p_intid, starts, interface_counts, grid.p_int.data()),
               "write variable 'p_int'");
  check_netcdf(nc_put_vara_double(ncid, p_midid, starts, counts, grid.p_mid.data()),
               "write variable 'p_mid'");
}

void NetcdfWriter::write_states(const std::vector<StateConst>& arrays, const std::size_t case_idx) {
  require_nonempty(arrays.size(), "state array list");
  const int caseid = get_dimension_id(ncid, "case");
  const int levid = get_dimension_id(ncid, "lev");
  int timeid;
  const int time_status = nc_inq_dimid(ncid, "time", &timeid);
  if (time_status == NC_EBADDIM) {
    check_netcdf(nc_def_dim(ncid, "time", arrays.size(), &timeid), "define dimension 'time'");
  } else {
    check_netcdf(time_status, "query dimension 'time'");
  }

  const int scalar_dimids[2] = {caseid, timeid};
  const int level_dimids[3] = {caseid, timeid, levid};
  const auto& var_descs = arrays.front().var_descs();
  std::vector<int> varids;
  varids.reserve(var_descs.size());
  for (const spaecies::VarDescPtr& var_desc : var_descs) {
    const bool is_scalar = var_desc->dimensions.empty();
    varids.push_back(get_or_define_variable(ncid, var_desc->name.c_str(), NC_DOUBLE,
                                            is_scalar ? 2 : 3,
                                            is_scalar ? scalar_dimids : level_dimids));
  }

  for (std::size_t time_idx = 0; time_idx != arrays.size(); ++time_idx) {
    for (std::size_t var_idx = 0; var_idx != varids.size(); ++var_idx) {
      auto var = arrays[time_idx].get_variable(var_descs[var_idx]->name).value();
      if (var_descs[var_idx]->dimensions.empty()) {
        const std::size_t starts[2] = {case_idx, time_idx};
        const std::size_t counts[2] = {1, 1};
        check_netcdf(nc_put_vara_double(ncid, varids[var_idx], starts, counts, &var[0]),
                     "write variable '" + var_descs[var_idx]->name + "'");
      } else {
        const std::size_t starts[3] = {case_idx, time_idx, 0};
        const std::size_t counts[3] = {1, 1, var.size()};
        check_netcdf(nc_put_vara_double(ncid, varids[var_idx], starts, counts, &var[0]),
                     "write variable '" + var_descs[var_idx]->name + "'");
      }
    }
  }
}

void NetcdfWriter::write_derived_vars(const std::vector<RainshaftDerivedVars>& dvars,
                                      const std::size_t case_idx) {
  require_nonempty(dvars.size(), "derived-variable list");
  const int caseid = get_dimension_id(ncid, "case");
  const int levid = get_dimension_id(ncid, "lev");
  const int ilevid = get_dimension_id(ncid, "ilev");
  const int timeid = get_dimension_id(ncid, "time");
  const int dimids[3] = {caseid, timeid, levid};
  const int interface_dimids[3] = {caseid, timeid, ilevid};
  const int z_intid = get_or_define_variable(ncid, "z_int", NC_DOUBLE, 3, interface_dimids);
  const int dzid = get_or_define_variable(ncid, "dz", NC_DOUBLE, 3, dimids);
  const int rho_dryid = get_or_define_variable(ncid, "rho_dry", NC_DOUBLE, 3, dimids);
  const int lambdarid = get_or_define_variable(ncid, "lambdar", NC_DOUBLE, 3, dimids);

  for (std::size_t time_idx = 0; time_idx != dvars.size(); ++time_idx) {
    const std::size_t nlev = dvars[time_idx].dz.size();
    const std::size_t starts[3] = {case_idx, time_idx, 0};
    const std::size_t counts[3] = {1, 1, nlev};
    const std::size_t interface_counts[3] = {1, 1, nlev + 1};
    check_netcdf(nc_put_vara_double(ncid, z_intid, starts, interface_counts, dvars[time_idx].z_int.data()),
                 "write variable 'z_int'");
    check_netcdf(nc_put_vara_double(ncid, dzid, starts, counts, dvars[time_idx].dz.data()),
                 "write variable 'dz'");
    check_netcdf(nc_put_vara_double(ncid, rho_dryid, starts, counts, dvars[time_idx].rho_dry.data()),
                 "write variable 'rho_dry'");
    check_netcdf(nc_put_vara_double(ncid, lambdarid, starts, counts, dvars[time_idx].lambdar.data()),
                 "write variable 'lambdar'");
  }
}

void NetcdfWriter::write_num_rhs_evals(const std::int64_t num_rhs_evals, const std::size_t case_idx) {
  const int caseid = get_dimension_id(ncid, "case");
  const int evalsid = get_or_define_variable(ncid, "num_rhs_evals", NC_INT64, 1, &caseid);
  check_netcdf(nc_put_var1(ncid, evalsid, &case_idx, &num_rhs_evals),
               "write variable 'num_rhs_evals'");
}

void NetcdfWriter::write_walltime_ms(const double walltime_ms, const std::size_t case_idx) {
  const int caseid = get_dimension_id(ncid, "case");
  const int walltime_msid = get_or_define_variable(ncid, "walltime_ms", NC_DOUBLE, 1, &caseid);
  check_netcdf(nc_put_var1_double(ncid, walltime_msid, &case_idx, &walltime_ms),
               "write variable 'walltime_ms'");
}

void NetcdfWriter::write_metadata(const int order, const double dt, const double dt_partition_1,
                                  const double dt_partition_2, const double rel_tol,
                                  const bool postprocess, const bool use_lookup,
                                  const std::string method_type, const int steps_per_output,
                                  const std::string initial_condition_file, const int num_cases,
                                  const int icase_in, const double final_time,
                                  const std::string processes) {
  const int orderid = define_scalar_variable(ncid, "method_order", NC_INT);
  check_netcdf(nc_put_var_int(ncid, orderid, &order), "write variable 'method_order'");

  const int dtid = define_scalar_variable(ncid, "dt", NC_DOUBLE);
  check_netcdf(nc_put_var_double(ncid, dtid, &dt), "write variable 'dt'");
  const int dt_partition_1id = define_scalar_variable(ncid, "dt_partition_1", NC_DOUBLE);
  check_netcdf(nc_put_var_double(ncid, dt_partition_1id, &dt_partition_1),
               "write variable 'dt_partition_1'");
  const int dt_partition_2id = define_scalar_variable(ncid, "dt_partition_2", NC_DOUBLE);
  check_netcdf(nc_put_var_double(ncid, dt_partition_2id, &dt_partition_2),
               "write variable 'dt_partition_2'");
  const int rel_tolid = define_scalar_variable(ncid, "rel_tol", NC_DOUBLE);
  check_netcdf(nc_put_var_double(ncid, rel_tolid, &rel_tol), "write variable 'rel_tol'");

  const auto postprocess_value = static_cast<unsigned char>(postprocess);
  const int postprocessid = define_scalar_variable(ncid, "postprocess", NC_UBYTE);
  check_netcdf(nc_put_var_uchar(ncid, postprocessid, &postprocess_value),
               "write variable 'postprocess'");
  const auto use_lookup_value = static_cast<unsigned char>(use_lookup);
  const int use_lookupid = define_scalar_variable(ncid, "use_lookup", NC_UBYTE);
  check_netcdf(nc_put_var_uchar(ncid, use_lookupid, &use_lookup_value),
               "write variable 'use_lookup'");

  const int method_typeid = define_scalar_variable(ncid, "method_type", NC_STRING);
  const char* method_type_value = method_type.c_str();
  check_netcdf(nc_put_var_string(ncid, method_typeid, &method_type_value),
               "write variable 'method_type'");

  const int steps_per_outputid = define_scalar_variable(ncid, "steps_per_output", NC_INT);
  check_netcdf(nc_put_var_int(ncid, steps_per_outputid, &steps_per_output),
               "write variable 'steps_per_output'");

  const int initial_condition_fileid = define_scalar_variable(ncid, "initial_condition_file", NC_STRING);
  const char* initial_condition_file_value = initial_condition_file.c_str();
  check_netcdf(nc_put_var_string(ncid, initial_condition_fileid, &initial_condition_file_value),
               "write variable 'initial_condition_file'");

  const int num_casesid = define_scalar_variable(ncid, "num_cases", NC_INT);
  check_netcdf(nc_put_var_int(ncid, num_casesid, &num_cases), "write variable 'num_cases'");
  const int icase_inid = define_scalar_variable(ncid, "icase_in", NC_INT);
  check_netcdf(nc_put_var_int(ncid, icase_inid, &icase_in), "write variable 'icase_in'");
  const int final_timeid = define_scalar_variable(ncid, "final_time", NC_DOUBLE);
  check_netcdf(nc_put_var_double(ncid, final_timeid, &final_time), "write variable 'final_time'");

  const int processesid = define_scalar_variable(ncid, "processes", NC_STRING);
  const char* processes_value = processes.c_str();
  check_netcdf(nc_put_var_string(ncid, processesid, &processes_value),
               "write variable 'processes'");
}
