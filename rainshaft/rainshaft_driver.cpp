#include "spaecies.hpp"
#include "evaporation.hpp"
#include "explicit_integrator.hpp"
#include "fixed_substep_integrator.hpp"
#include "nudging.hpp"
#include "operator_splitting_integrator.hpp"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_solution.hpp"
#include "rainshaft_tend_descs.hpp"
#include "rainshaft_ncio.hpp"
#include "rainshaft_sum_process.hpp"
#include "saturation.hpp"
#include "sed_cfl_integrator.hpp"
#include "sedimentation.hpp"
#include "self_collision.hpp"
#include "sequential_split_integrator.hpp"
#include <iostream>
#include <cmath>
#include <chrono>
#include <memory>
#include <string>
#include "imex_integrator.hpp"
#include "mri_integrator.hpp"
#include <boost/program_options.hpp>

int main(int argc, char* argv[])
{
  using std::chrono::high_resolution_clock;
  using std::chrono::duration;
  namespace po = boost::program_options;

  // Set up model constants.
  // SPS: Choose rho_top in a more principled way?
  RainshaftConstants constants{3.14159265358979323846,
                               287.04, 1.00464e3, 461.50, 997., 2.501e6,
                               0.62197, 1.e-14, 9.80616, 1.e-5, 5.e-3,
                               0.988919555598356, 1.e3, 1.e-4};
  // Parse input
  int order;
  double dt;
  double dt_slow, dt_slow_in;
  std::string name;
  int steps_per_output;
  std::string simulation_name;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("order", po::value<int>(&order)->default_value(2), "order of method")
		("dt", po::value<double>(&dt)->default_value(0.1), "step size. if integration type is set to MRI, this argument is the inner time step size")
    ("dt_slow", po::value<double>(&dt_slow_in), "slow/outer step size for MRI integration. negative values indicate ratios, e.g. dt_slow=-0.5 corresponds to dt_slow = dt_fast/2")
		("type", po::value<std::string>(&name)->default_value("explicit"), "type of integrator (e.g. explicit, implicit, imex, mri)")
    ("steps", po::value<int>(&steps_per_output)->default_value(-1), "frequency of saved output, e.g. write output to netCDF every steps_per_output timesteps. -1 to save only the first and last steps")
    ("simname", po::value<std::string>(), "common name for group of simulations")
  ;

	po::variables_map vm{};
	po::store(parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
		std::cout << desc << "\n";
		return 1;
	}
 
	if (vm["order"].defaulted()) {
    std::cout << "Order of integration method was not set. Defaulting to " << order << "." << std::endl;
  } else {
    std::cout << "Order of integration method set to " << order << "." << std::endl;
  }

	if (vm["dt"].defaulted()) {
    std::cout << "Time step was not set. Defaulting to " << dt << "." << std::endl;
  } else {
    std::cout << "Time step was set to " << dt << "." << std::endl;
  }

  if (vm.count("dt_slow")) {
    dt_slow_in = vm["dt_slow"].as<double>();
    if (dt_slow_in < 0.0) {
      dt_slow = abs(dt_slow_in) * dt;
    } else {
      dt_slow = dt_slow_in;
    }
    std::cout << "Slow/outer time step was set to " << dt_slow << "." << std::endl;
  } else {
    std::cout << "Slow/outer time step was not set." << std::endl;
  }

	if (vm["type"].defaulted()) {
    std::cout << "Integration type was not set. Defaulting to " << name << "." << std::endl;
  } else {
    std::cout << "Integration type was set to " << name << "." << std::endl;
  }

  if (vm["steps"].defaulted()) {
    std::cout << "Save frequency not set. Defaulting to " << steps_per_output << "." << std::endl;
  } else {
    std::cout << "Save frequency set to " << steps_per_output << "." << std::endl;
  }

  if (vm.count("simname")) {
    simulation_name = vm["simname"].as<std::string>();
    std::cout << "Simulation name: " << simulation_name << "." << std::endl;
  }
  
  // Approximate model top in meters.
  // (The grid maker will actually use the next higher-altitude E3SM level.)
  double model_top = 2.e3;
  // Surface pressure and temperature in Pa and K, respectively.
  double srf_pres = 1.e5, srf_temp = 293.15;
  // Lapse rate of initial condition in K/m.
  double lapse_rate = 6.5e-3;
  // Relative humidity of initial condition.
  double rel_hum_init = 0.7;
  // Time scale over which to nudge t and q back to initial condition in seconds.
  double nudge_time_scale = 15. * 60.;
  // Time of simulation start.
  double initial_time = 0.;
  // Final time to integrate to.
  double final_time = 1800.;
  RainshaftGrid grid = make_e3sm_like_grid(constants, model_top, srf_pres,
                                           srf_temp, lapse_rate);
  std::size_t nlev = grid.nlev;
  // Set up initial condition.
  SaturationFormulae sat_form(constants);

  spaecies::Domain dom;
  spaecies::DimensionPtr lev_dim = dom.add_dimension("level", nlev);
  spaecies::VarDescPtr t_desc = dom.add_var_desc("T", spaecies::Float64Type, {lev_dim}, "K");
  spaecies::VarDescPtr q_desc = dom.add_var_desc("q", spaecies::Float64Type, {lev_dim}, "kg/kg");
  spaecies::VarDescPtr nr_desc = dom.add_var_desc("nr", spaecies::Float64Type, {lev_dim}, "1/kg");
  spaecies::VarDescPtr qr_desc = dom.add_var_desc("qr", spaecies::Float64Type, {lev_dim}, "kg/kg");
  VarDescList state_descs = {t_desc, q_desc, nr_desc, qr_desc};
  VarDescList tend_descs = tend_descs_from_state_descs(dom, state_descs);
  State initial_state(state_descs);
  VarMut t = initial_state.get_variable("T");
  VarMut q = initial_state.get_variable("q");
  VarMut nr = initial_state.get_variable("nr");
  VarMut qr = initial_state.get_variable("qr");
  for (std::size_t i = 0; i != nlev; ++i) {
    nr[i] = 0.;
    qr[i] = 0.;
  }

  // Coming up with an initial condition for t and q is slightly tricky, because
  // we only have an implicit relationship between t, q, and dz. But since the
  // effect of q on layer height is not large, start by ignoring it, in which
  // case we do have an explicit relationship between t and dz.
  double rog = constants.rdry / constants.g;
  std::vector<double> z_int(nlev+1, 0.);
  for (int il = nlev - 1; il != -1; --il) {
    double pdel = grid.p_int[il+1]-grid.p_int[il];
    t[il] = (srf_temp - lapse_rate*z_int[il+1])
      / (1. + 0.5 * lapse_rate * rog * pdel/ grid.p_mid[il]);
    z_int[il] = rog * t[il] * pdel / grid.p_mid[il];
  }
  // Now iterate until change in dz between successive iterations is small;
  // ten iterations should do it.
  bool converged = false;
  double t_tolerance = 1.e-3; // Kelvin
  for (int it = 0; it != 10 && !converged; ++it) {
    std::vector<double> t_v(nlev, 0.);
    for (std::size_t il = 0; il != nlev; ++il) {
      q[il] = rel_hum_init * sat_form.q_sat_dry(t[il], grid.p_mid[il]);
      // Virtual temperature factor.
      double t_v_fac = 1. + ((1/constants.epsilon_h2o - 1.) * q[il]);
      t_v[il] = t[il] * t_v_fac;
    }
    auto dz = grid.calc_dz(constants, t_v);
    z_int = dz_to_z_int(dz);
    converged = true;
    for (std::size_t il = 0; il!= nlev; ++il) {
      double new_temp = srf_temp - lapse_rate * (z_int[il] - 0.5 * dz[il]);
      if (abs(new_temp - t[il]) > t_tolerance) {
        converged = false;
      }
      t[il] = new_temp;
    }
  }
  if (!converged) {
    std::cerr << "Initial condition iteration failed to converge." << std::endl;
    return 1;
  }
  RainshaftDerivedVars initial_dvars = RainshaftDerivedVars(constants, grid, initial_state);
  // Sedimentation process.
  Sedimentation sed(constants, true, false);
  // Self-collision processes.
  SelfCollision self_coll;
  // Evaporation process.
  Evaporation evap(constants, sat_form, true, false);
  // Nudging to initial condition.
  // SPS: Need some kind of span-like interface to avoid having to
  // do this copy.
  std::vector<double> t_vec, q_vec;
  for (std::size_t i = 0; i != nlev; ++i) {
    t_vec.push_back(t[i]);
    q_vec.push_back(q[i]);
  }
  Nudging nudge(nudge_time_scale, t_vec, q_vec);
  // Sum of all processes.
  SumProcess exp_processes = SumProcess{{&evap, &nudge, &self_coll}};
  // Sum of local processes.
  SumProcess imp_processes = SumProcess{{&sed}};

  SumProcess all_processes = SumProcess{{&sed, &nudge, &self_coll, &evap}};

  const auto intg = [&]() -> std::unique_ptr<RainshaftIntegrator> {
    if (name == "explicit") {
      return std::make_unique<ExplicitIntegrator>(constants, grid, &all_processes, state_descs, tend_descs, dt, order, steps_per_output);
    } else if (name == "implicit") {
      return std::make_unique<IMEXIntegrator>(constants, grid, nullptr, &all_processes, state_descs, tend_descs, dt, order, steps_per_output);
    } else if (name == "imex") {
      return std::make_unique<IMEXIntegrator>(constants, grid, &exp_processes, &imp_processes, state_descs, tend_descs, dt, order, steps_per_output);
    } else if (name == "mri") {
      return std::make_unique<MRIIntegrator>(constants, grid, &imp_processes, &exp_processes, nullptr, state_descs, tend_descs, dt, dt_slow, order, steps_per_output);
    } else if (name == "split") {
      return std::make_unique<OperatorSplittingIntegrator>(constants, grid, &exp_processes, &imp_processes, state_descs, tend_descs, dt_slow, dt, dt, order, steps_per_output);
    } else {
      throw std::logic_error("Invalid integration method. Options: explicit, imex, mri.");
    }
  }();

  auto before_sol = high_resolution_clock::now();
  RainshaftSolution solution = intg->integrate(initial_time, final_time, initial_state);
  auto after_sol = high_resolution_clock::now();
  // Time taken for solution.
  duration<double, std::milli> walltime_ms = after_sol - before_sol;
  std::cout << "Time: " << walltime_ms.count() << std::endl;

  // parse time step for savefile name
  std::string str1 = std::to_string(dt);
  str1.erase(str1.find_first_of('.'), 1);
  str1.erase(0, str1.find_first_not_of('0'));
  str1.erase ( str1.find_last_not_of('0') + 1, std::string::npos );

  std::string filename;
  // Write out grid and all states.
  if (vm.count("simname")) {
    if (vm.count("dt_slow")) {
      filename = "results/rainshaft_" + vm["simname"].as<std::string>() + "_order" + std::to_string(order) + "_dtslow" + std::to_string(dt_slow_in) + "_dt" + std::to_string(dt) + ".nc";
    } else {
      filename = "results/rainshaft_" + vm["simname"].as<std::string>() + "_order" + std::to_string(order) + "_dt" + std::to_string(dt) + ".nc";
    }
  } else {
    if (vm.count("dt_slow")) {
      filename = "results/rainshaft_order" + std::to_string(order) + "_dtslow" + std::to_string(dt_slow_in) +  "_dt" + std::to_string(dt) + ".nc";
    } else {
      filename = "results/rainshaft_order" + std::to_string(order) +  "_dt" + std::to_string(dt) + ".nc";
    }
  }

  NetcdfWriter writer(filename);
  writer.write_grid(grid);
  writer.write_states(solution.states);
  std::vector<RainshaftDerivedVars> solution_dvars;
  solution_dvars.reserve(solution.states.size());
  for (StateConst state : solution.states) {
    solution_dvars.emplace_back(constants, grid, state);
  }
  writer.write_derived_vars(solution_dvars);
  writer.write_num_rhs_evals(solution.num_rhs_evals);
  writer.write_walltime_ms(walltime_ms.count());
  // Ensure that the library is linked and greet the user.
  spaecies::do_nothing();
  return 0;
}
