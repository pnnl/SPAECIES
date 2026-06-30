#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <boost/program_options.hpp>

#include "libspaecies/spaecies.hpp"
#include "rainshaft/rainshaft_types.hpp"

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
  std::vector<double> read_grid(std::size_t case_idx, std::size_t time) const;
  void read_state(std::size_t case_idx, State &initial_state, std::size_t time) const;
};


int main(int argc, char* argv[])
{
  namespace po = boost::program_options;
  using VarMut = spaecies::ContiguousVariableView<double>;

  // Command line arguments
  std::string input_file;
  std::size_t time;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("input,i", po::value(&input_file)->required(), "input netCDF file")
    ("time,t", po::value(&time)->default_value(0), "time index")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  NetcdfReader reader(input_file);
  std::size_t num_cases, num_levs;
  std::tie(num_cases, num_levs) = reader.read_num_cases_and_levs();

  // TODO: consider creating a helper function in rainshaft_driver that creates a RainshaftState
  spaecies::Domain dom;
  spaecies::DimensionPtr lev_dim = dom.add_dimension("level", num_levs);
  spaecies::VarDescPtr T_desc = dom.add_var_desc("T", spaecies::Float64Type, {lev_dim}, "K");
  spaecies::VarDescPtr q_desc = dom.add_var_desc("q", spaecies::Float64Type, {lev_dim}, "kg/kg");
  spaecies::VarDescPtr nr_desc = dom.add_var_desc("nr", spaecies::Float64Type, {lev_dim}, "1/kg");
  spaecies::VarDescPtr qr_desc = dom.add_var_desc("qr", spaecies::Float64Type, {lev_dim}, "kg/kg");
  State state({T_desc, q_desc, nr_desc, qr_desc});

  for (std::size_t icase=0; icase < num_cases; icase++)
  {
    std::vector<double> z_int = reader.read_grid(icase, time);
    reader.read_state(icase, state, time);
    std::cout << state.get_variable("T")[0] << std::endl;
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

std::vector<double> NetcdfReader::read_grid(std::size_t case_idx, std::size_t time) const
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

void NetcdfReader::read_state(std::size_t case_idx, State &initial_state, std::size_t time) const
{
  int lev_id, T_id, q_id, nr_id, qr_id;
  nc_inq_varid(ncid, "lev", &lev_id);
  nc_inq_varid(ncid, "T", &T_id);
  nc_inq_varid(ncid, "q", &q_id);
  nc_inq_varid(ncid, "nr", &nr_id);
  nc_inq_varid(ncid, "qr", &qr_id);
  std::size_t num_levs;
  nc_inq_dimlen(ncid, lev_id, &num_levs);
  std::size_t starts[3] = {case_idx, time, 0};
  std::size_t counts[3] = {1, 1, num_levs};
  spaecies::ContiguousVariableView<double> T = initial_state.get_variable("T");
  spaecies::ContiguousVariableView<double> q = initial_state.get_variable("q");
  spaecies::ContiguousVariableView<double> nr = initial_state.get_variable("nr");
  spaecies::ContiguousVariableView<double> qr = initial_state.get_variable("qr");
  nc_get_vara_double(ncid, T_id, starts, counts, &T[0]);
  nc_get_vara_double(ncid, q_id, starts, counts, &q[0]);
  nc_get_vara_double(ncid, nr_id, starts, counts, &nr[0]);
  nc_get_vara_double(ncid, qr_id, starts, counts, &qr[0]);
}
