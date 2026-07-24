#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <optional>
#include <type_traits>
#include <boost/program_options.hpp>
#include <matplot/matplot.h>

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

class NetcdfFile
{
  int file_id;
  const std::string file_name;
public:
  NetcdfFile(const std::string& file_name);
  ~NetcdfFile();
  NetcdfFile(const NetcdfFile&) = delete;
  NetcdfFile(NetcdfFile&&) = delete;
  NetcdfFile& operator=(const NetcdfFile&) = delete;
  NetcdfFile& operator=(NetcdfFile&&) = delete;
  // CJV: revisit whether this should be const given the underlying file stream
  //      is likely changed (even if the overall operation is read-only)
  int get_file_id() const { return file_id; }
  int get_dim_id(const std::string& dim_name) const;
  int get_var_id(const std::string& var_name) const;
  template <typename T>
  T get_var_val(const std::string& var_name, const std::size_t idx) const;
};

// variant NetcdfReader class modified for output of NetcdfWriter class
class NetcdfReader
{
  NetcdfFile file;
public:
  NetcdfReader(const std::string& file_name) : file(file_name) {}
  std::size_t read_num_cases() const;
  RainshaftGrid read_grid_pressure(const std::size_t case_idx) const;
  std::vector<double> read_grid_height(const std::size_t case_idx, const std::size_t time_idx) const;
  void read_boundary_conditions(const std::size_t case_idx, RainshaftConstants &constants) const;
  void read_state(const std::size_t case_idx, const std::size_t time_idx, State &initial_state) const;
};

double calc_max_characteristic_speed(const RainshaftConstants& constants, const double lambdar);


int main(int argc, char* argv[])
{
  const double epsilon = std::numeric_limits<double>::epsilon();

  // parse input file name and time index from command line arguments
  namespace po = boost::program_options;
  std::string input_file;
  std::size_t time_idx;
  bool use_zero_mur;
  int case_idx;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("input,i", po::value(&input_file)->required(), "input netCDF file")
    ("time_idx,ti", po::value(&time_idx)->default_value(0), "time index")
    ("use_zero_mur", po::value(&use_zero_mur)->default_value(defaults::use_zero_mur), "use zero for rain shape parameter mu (legacy value)")
    ("case_idx", po::value(&case_idx)->default_value(-1), "(optional) specific case to analyze (otherwise all cases are analyzed).")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help"))
	{
		std::cout << desc << "\n";
		return 1;
	}
  po::notify(vm);

  const NetcdfReader reader(input_file);
  std::size_t num_cases = reader.read_num_cases();

  // create rainshaft processes
  const double mur = use_zero_mur ? 0.0 : 1.0;
  RainshaftConstants constants = create_RainshaftConstants(defaults::qsmall,
    defaults::epsilon_qsat_fac, defaults::epsilon_self_coll, mur);
  SaturationFormulae sat_form(constants);
  Evaporation evap(constants, sat_form, defaults::use_lookup, false,
                   defaults::regularize_qsat, std::nullopt);
  SelfCollision self_coll(defaults::regularize_lambdar);
  Sedimentation sed(constants, defaults::use_lookup, false);

  // compute process rates for each case in input file
  std::size_t case_start = 0;
  if (case_idx >= 0)
  {
    case_start = (std::size_t) case_idx;
    num_cases = 1;
  }
  std::vector<double> evaporation_rates(num_cases);
  std::vector<double> self_coll_rates(num_cases);
  std::vector<double> sedimentation_accuracy_rates(num_cases);
  std::vector<double> sedimentation_stability_rates(num_cases);
  for (std::size_t icase=case_start; icase < case_start + num_cases; icase++)
  {
    // read in grid and boundary condition information for current case
    reader.read_boundary_conditions(icase, constants);
    const std::vector<double> z = reader.read_grid_height(icase, time_idx);
    RainshaftGrid grid = reader.read_grid_pressure(icase);

    // create state and tendency objects
    spaecies::Domain dom;
    spaecies::DimensionPtr lev_dim = dom.add_dimension("level", grid.nlev);
    VarDescList state_descs = get_prognostic_variables(dom, lev_dim);
    State state(state_descs);
    VarDescList tend_descs = tend_descs_from_state_descs(dom, state_descs);
    Tendency evap_tend(tend_descs), self_coll_tend(tend_descs), sed_tend(tend_descs);

    // read in state and compute derived variables and tendencies
    reader.read_state(icase, time_idx, state);
    RainshaftDerivedVars dvars(constants, grid, state, defaults::regularize_lambdar);
    evap.calc_tend(constants, grid, state, dvars, evap_tend);
    self_coll.calc_tend(constants, grid, state, dvars, self_coll_tend);
    sed.calc_tend(constants, grid, state, dvars, sed_tend);

    // create helper views
    spaecies::ContiguousVariableView<const double> T = state.get_variable("T").value();
    spaecies::ContiguousVariableView<const double> q = state.get_variable("q").value();
    spaecies::ContiguousVariableView<const double> nr = state.get_variable("nr").value();
    spaecies::ContiguousVariableView<const double> qr = state.get_variable("qr").value();
    spaecies::ContiguousVariableView<const double> Sqr = evap_tend.get_variable("qr_tend").value();
    spaecies::ContiguousVariableView<const double> Srsc = self_coll_tend.get_variable("nr_tend").value();

    // compute rates
    double evaporation_rate = 0.0;
    double self_coll_rate = 0.0;
    double sedimentation_accuracy_rate = 0.0;
    double sedimentation_stability_rate = 0.0;
    for (std::size_t ilev = 0; ilev != grid.nlev; ilev++)
    {
      if (qr[ilev] > epsilon)
      {
        evaporation_rate = std::max(evaporation_rate, std::abs(Sqr[ilev]) / qr[ilev]);
      }
      if (nr[ilev] > epsilon)
      {
        self_coll_rate = std::max(self_coll_rate, std::abs(Srsc[ilev]) / nr[ilev]);
      }
      const double sigma = calc_max_characteristic_speed(constants, dvars.lambdar[ilev]);
      if (0 < ilev && ilev < grid.nlev && qr[ilev] > epsilon)
      {
        sedimentation_accuracy_rate = std::max(sedimentation_accuracy_rate, sigma /
          std::sqrt( qr[ilev] /
            (std::abs(((qr[ilev-1] - qr[ilev]) / (z[ilev-1] - z[ilev]) -
                       (qr[ilev] - qr[ilev+1]) / (z[ilev] - z[ilev+1])) /
                      (0.5*(z[ilev-1] - z[ilev+1]))) + epsilon)
                   ));
      }
      sedimentation_stability_rate = std::max(sedimentation_stability_rate,
        sigma / (z[ilev] - z[ilev+1]));
    }
    const std::size_t index = icase - case_start;
    evaporation_rates[index] = evaporation_rate;
    self_coll_rates[index] = self_coll_rate;
    sedimentation_accuracy_rates[index] = sedimentation_accuracy_rate;
    sedimentation_stability_rates[index] = sedimentation_stability_rate;
  }

  matplot::line_handle plot;
  matplot::hold(matplot::on);
  //matplot::line_handle plot = matplot::plot(evaporation_rates, "o");
  //plot->display_name("evaporation rate");
  //plot = matplot::plot(self_coll_rates, "o");
  //plot->display_name("self collection rate");
  plot = matplot::plot(sedimentation_accuracy_rates, "o");
  plot->display_name("sedimentation accuracy rate");
  plot = matplot::plot(sedimentation_stability_rates, "o");
  plot->display_name("sedimentation stability rate");
  matplot::hold(matplot::off);

  matplot::xlabel("case number");
  matplot::ylabel("process rate / inverse timescale (1/s)");
  matplot::legend();

  matplot::show();

  return 0;
}

NetcdfFile::NetcdfFile(const std::string& file_name) : file_name(file_name)
{
  // CJV: check if there should be an error check here
  nc_open(file_name.c_str(), NC_NOWRITE, &file_id);
}

NetcdfFile::~NetcdfFile()
{
  nc_close(file_id);
}

int NetcdfFile::get_dim_id(const std::string& dim_name) const
{
  int dim_id;
  if (nc_inq_dimid(file_id, dim_name.c_str(), &dim_id) == NC_ENOTVAR)
  {
    throw std::runtime_error("NetCDF file " + file_name + " does not contain dimension " + dim_name);
  }
  return dim_id;
}

int NetcdfFile::get_var_id(const std::string& var_name) const
{
  int var_id;
  if (nc_inq_varid(file_id, var_name.c_str(), &var_id) == NC_ENOTVAR)
  {
    throw std::runtime_error("NetCDF file " + file_name + " does not contain variable " + var_name);
  }
  return var_id;
}

template <typename T>
T NetcdfFile::get_var_val(const std::string& var_name, const std::size_t idx) const
{
  T value;
  const int var_id = get_var_id(var_name);
  if constexpr (std::is_integral_v<T>)
    nc_get_var1_int(file_id, var_id, &idx, &value);
  else
    static_assert(false, "Unsupported type");

  return value;
}

std::size_t NetcdfReader::read_num_cases() const
{
  const int case_id = file.get_dim_id("case");
  std::size_t num_cases;
  nc_inq_dimlen(file.get_file_id(), case_id, &num_cases);
  return num_cases;
}

std::vector<double> NetcdfReader::read_grid_height(const std::size_t case_idx, const std::size_t time_idx) const
{
  const int num_ilevs = file.get_var_val<int>("nlev", case_idx) + 1;
  const int z_int_id = file.get_var_id("z_int");
  std::size_t starts[3] = {case_idx, time_idx, 0};
  std::size_t counts[3] = {1, 1, (std::size_t) num_ilevs};
  std::vector<double> z_int(num_ilevs);
  nc_get_vara_double(file.get_file_id(), z_int_id, starts, counts, z_int.data());
  return z_int;
}

RainshaftGrid NetcdfReader::read_grid_pressure(const std::size_t case_idx) const
{
  const int num_ilevs = file.get_var_val<int>("nlev", case_idx) + 1;
  const int p_int_id = file.get_var_id("p_int");
  const std::size_t starts[2] = {case_idx, 0};
  const std::size_t counts[2] = {1, (std::size_t) (num_ilevs)};
  std::vector<double> p_int(num_ilevs);
  nc_get_vara_double(file.get_file_id(), p_int_id, starts, counts, p_int.data());
  return RainshaftGrid(p_int);
}

void NetcdfReader::read_boundary_conditions(const std::size_t case_idx, RainshaftConstants &constants) const
{
  const int rho_top_id = file.get_var_id("rho_top");
  const int qr_top_id = file.get_var_id("qr_top");
  const int nr_top_id = file.get_var_id("nr_top");
  nc_get_var1_double(file.get_file_id(), rho_top_id, &case_idx, &constants.rho_top);
  nc_get_var1_double(file.get_file_id(), qr_top_id, &case_idx, &constants.qr_top);
  nc_get_var1_double(file.get_file_id(), nr_top_id, &case_idx, &constants.nr_top);
}

void NetcdfReader::read_state(const std::size_t case_idx, const std::size_t time_idx, State &state) const
{
  const int num_levs = file.get_var_val<int>("nlev", case_idx);
  const int T_id = file.get_var_id("T");
  const int q_id = file.get_var_id("q");
  const int nr_id = file.get_var_id("nr");
  const int qr_id = file.get_var_id("qr");
  std::size_t starts[3] = {case_idx, time_idx, 0};
  std::size_t counts[3] = {1, 1, (std::size_t) num_levs};
  spaecies::ContiguousVariableView<double> T = state.get_variable("T").value();
  spaecies::ContiguousVariableView<double> q = state.get_variable("q").value();
  spaecies::ContiguousVariableView<double> nr = state.get_variable("nr").value();
  spaecies::ContiguousVariableView<double> qr = state.get_variable("qr").value();
  nc_get_vara_double(file.get_file_id(), T_id, starts, counts, &T[0]);
  nc_get_vara_double(file.get_file_id(), q_id, starts, counts, &q[0]);
  nc_get_vara_double(file.get_file_id(), nr_id, starts, counts, &nr[0]);
  nc_get_vara_double(file.get_file_id(), qr_id, starts, counts, &qr[0]);
}

double calc_max_characteristic_speed(const RainshaftConstants& constants, const double lambdar)
{
  const auto [vnr_info, vqr_info] =
    Sedimentation::rain_fall_speeds_gamma<true>(constants, lambdar);
  const double vnr = get_val(vnr_info);
  const double dvnr_dlam = get_grad(vnr_info)[0];
  const double vqr = get_val(vqr_info);
  const double dvqr_dlam = get_grad(vqr_info)[0];
  const double vnr_term = vnr + dvnr_dlam * lambdar/3.0;
  const double vqr_term = vqr - dvqr_dlam * lambdar/3.0;
  const double b = -(vnr_term + vqr_term);
  const double c = vnr_term*vqr_term + dvnr_dlam*dvqr_dlam*std::pow(lambdar/3.0,2);
  return 0.5*(-b + std::sqrt(b*b - 4.0*c));
}
