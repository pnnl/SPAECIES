#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <optional>
#include <boost/program_options.hpp>

#include "libspaecies/spaecies.hpp"
#include "rainshaft/evaporation.hpp"
#include "rainshaft/rainshaft_driver.hpp"
#include "rainshaft/rainshaft_types.hpp"
#include "rainshaft/saturation.hpp"
#include "rainshaft/sedimentation.hpp"
#include "rainshaft/self_collision.hpp"
#include "rainshaft/rainshaft_tend_descs.hpp"

extern "C"
{
#include "netcdf.h"
}


// variant NetcdfReader class modified for output of NetcdfWriter class
class NetcdfReader
{
  int ncid;
public:
  NetcdfReader(const std::string& file_name);
  ~NetcdfReader();
  NetcdfReader(const NetcdfReader&) = delete;
  NetcdfReader(NetcdfReader&&) = delete;
  NetcdfReader& operator=(const NetcdfReader&) = delete;
  NetcdfReader& operator=(NetcdfReader&&) = delete;
  std::tuple<std::size_t, std::size_t> read_num_cases_and_levs() const;
  std::vector<double> read_grid(std::size_t case_idx) const;
  void read_state(std::size_t case_idx, State &initial_state, std::size_t time_idx) const;
};


int main(int argc, char* argv[])
{
  const double epsilon = std::numeric_limits<double>::epsilon();

  // parse input file name and time index from command line arguments
  namespace po = boost::program_options;
  std::string input_file;
  std::size_t time_idx;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("input,i", po::value(&input_file)->required(), "input netCDF file")
    ("time_idx,ti", po::value(&time_idx)->default_value(0), "time index")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  const NetcdfReader reader(input_file);
  std::size_t num_cases, num_levs;
  std::tie(num_cases, num_levs) = reader.read_num_cases_and_levs();

  // create rainshaft processes
  RainshaftConstants constants = create_RainshaftConstants(defaults::qsmall,
    defaults::epsilon_qsat_fac, defaults::epsilon_self_coll);
  SaturationFormulae sat_form(constants);
  Evaporation evap(constants, sat_form, defaults::use_lookup, false,
                   defaults::regularize_qsat, std::nullopt);
  SelfCollision self_coll(defaults::regularize_lambdar);
  Sedimentation sed(constants, defaults::use_lookup, false);

  // create state and tendency objects
  spaecies::Domain dom;
  spaecies::DimensionPtr lev_dim = dom.add_dimension("level", num_levs);
  VarDescList state_descs = get_prognostic_variables(dom, lev_dim);
  State state(state_descs);
  VarDescList tend_descs = tend_descs_from_state_descs(dom, state_descs);
  Tendency evap_tend(tend_descs), self_coll_tend(tend_descs), sed_tend(tend_descs);

  // create helper views
  spaecies::ContiguousVariableView<const double> T = state.get_variable("T").value();
  spaecies::ContiguousVariableView<const double> q = state.get_variable("q").value();
  spaecies::ContiguousVariableView<const double> nr = state.get_variable("nr").value();
  spaecies::ContiguousVariableView<const double> qr = state.get_variable("qr").value();

  // compute process rates for each case in input file
  for (std::size_t icase=0; icase < num_cases; icase++)
  {

    reader.read_state(icase, state, time_idx);
    RainshaftGrid grid = reader.read_grid(icase);
    RainshaftDerivedVars dvars(constants, grid, state, defaults::regularize_lambdar);

    evap.calc_tend(constants, grid, state, dvars, evap_tend);
    self_coll.calc_tend(constants, grid, state, dvars, self_coll_tend);
    sed.calc_tend(constants, grid, state, dvars, sed_tend);

    double evaporation_rate = 0.0;
    double self_coll_rate = 0.0;
    double sedimentation_accurate_rate = 0.0;
    double sedimentation_stability_rate = 0.0;
    for (std::size_t ilev = 0; ilev != grid.nlev; ilev++)
    {
      evaporation_rate = std::max(evaporation_rate,
        qr[ilev] / (evap_tend.get_variable("qr_tend").value()[ilev] + epsilon));
      self_coll_rate = std::max(self_coll_rate,
        nr[ilev] / (self_coll_tend.get_variable("nr_tend").value()[ilev] + epsilon));
      }
    std::cout << "Evaporation: " << evaporation_rate << std::endl;
    std::cout << "Self-collection: " << self_coll_rate << std::endl;
  }

  return 0;
}

NetcdfReader::NetcdfReader(const std::string& file_name)
{
  nc_open(file_name.c_str(), NC_NOWRITE, &ncid);
}

NetcdfReader::~NetcdfReader()
{
  nc_close(ncid);
}

std::tuple<std::size_t, std::size_t> NetcdfReader::read_num_cases_and_levs() const
{
  int case_id, lev_id;
  nc_inq_dimid(ncid, "case", &case_id);
  nc_inq_dimid(ncid, "lev", &lev_id);
  std::size_t num_cases, num_levs;
  nc_inq_dimlen(ncid, case_id, &num_cases);
  nc_inq_dimlen(ncid, lev_id, &num_levs);
  return {num_cases, num_levs};
}

std::vector<double> NetcdfReader::read_grid(std::size_t case_idx) const
{
  int ilev_id, z_int_id;
  nc_inq_dimid(ncid, "ilev", &ilev_id);
  nc_inq_varid(ncid, "z_int", &z_int_id);
  std::size_t num_ilevs;
  nc_inq_dimlen(ncid, ilev_id, &num_ilevs);
  std::size_t starts[2] = {case_idx, 0};
  std::size_t counts[2] = {1, num_ilevs};
  std::vector<double> z_int(num_ilevs);
  nc_get_vara_double(ncid, z_int_id, starts, counts, z_int.data());
  return std::move(z_int);
}

void NetcdfReader::read_state(std::size_t case_idx, State &state, std::size_t time_idx) const
{
  int lev_id, T_id, q_id, nr_id, qr_id;
  nc_inq_varid(ncid, "lev", &lev_id);
  nc_inq_varid(ncid, "T", &T_id);
  nc_inq_varid(ncid, "q", &q_id);
  nc_inq_varid(ncid, "nr", &nr_id);
  nc_inq_varid(ncid, "qr", &qr_id);
  std::size_t num_levs;
  nc_inq_dimlen(ncid, lev_id, &num_levs);
  std::size_t starts[3] = {case_idx, time_idx, 0};
  std::size_t counts[3] = {1, 1, num_levs};
  spaecies::ContiguousVariableView<double> T = state.get_variable("T").value();
  spaecies::ContiguousVariableView<double> q = state.get_variable("q").value();
  spaecies::ContiguousVariableView<double> nr = state.get_variable("nr").value();
  spaecies::ContiguousVariableView<double> qr = state.get_variable("qr").value();
  nc_get_vara_double(ncid, T_id, starts, counts, &T[0]);
  nc_get_vara_double(ncid, q_id, starts, counts, &q[0]);
  nc_get_vara_double(ncid, nr_id, starts, counts, &nr[0]);
  nc_get_vara_double(ncid, qr_id, starts, counts, &qr[0]);
}
