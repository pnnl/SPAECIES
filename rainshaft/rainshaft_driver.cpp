#include "spaecies.hpp"
#include "evaporation.hpp"
#include "explicit_integrator.hpp"
#include "fixed_substep_integrator.hpp"
#include "nudging.hpp"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_initial.hpp"
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
#include <optional>
#include <string>
#include "imex_integrator.hpp"
#include "mri_integrator.hpp"

int main(int, char* argv[])
{
  using std::chrono::high_resolution_clock;
  using std::chrono::duration;

  // Get command line arguments
  std::size_t i = 0;
  std::string initial_condition(argv[++i]); // File for initial conditions, or 'adiabatic'.
  double sim_len = std::stod(argv[++i]); // Length of time to integrate to
  bool do_nudging;
  std::string nudge_opt(argv[++i]); // Whether or not to do nudging
  if (nudge_opt == "TRUE" || nudge_opt == "True" || nudge_opt == "T"
      || nudge_opt == "true" || nudge_opt == "t") {
    do_nudging = true;
  } else if (nudge_opt == "FALSE" || nudge_opt == "False" || nudge_opt == "F"
      || nudge_opt == "false" || nudge_opt == "f") {
    do_nudging = false;
  } else {
    throw std::logic_error("Unrecognized nudging option");
  }
  std::string method_type(argv[++i]); // Time integration method type
  std::optional<double> dt_fast; // Fast time step for MRI method
  if (method_type == "mri") {
    dt_fast = std::stod(argv[++i]);
  }
  double dt = std::stod(argv[++i]); // Overall time step size
  int order = std::stoi(argv[++i]); // Order of accuracy desired
  int steps_per_output = std::stoi(argv[++i]); // Number of time steps between each output
  std::string output_file(argv[++i]); // Name of file to output to

  // Set up model constants.
  // SPS: Choose rho_top in a more principled way?
  RainshaftConstants constants{3.14159265358979323846,
                               287.04, 1.00464e3, 461.50, 997., 2.501e6,
                               0.62197, 1.e-14, 9.80616, 1.e-5, 5.e-3,
                               0.988919555598356, 1.e3, 1.e-4};
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
  double final_time = initial_time + sim_len;
  // Saturation vapor pressure formula.
  SaturationFormulae sat_form(constants);

  // Set up max grid size
  std::size_t num_cases, max_levs;
  RainshaftGrid default_grid = make_e3sm_like_grid(constants, model_top, srf_pres,
                                                   srf_temp, lapse_rate);
  if (initial_condition == "adiabatic") {
    num_cases = 1;
    max_levs = default_grid.nlev;
  }
  NetcdfWriter writer(output_file, num_cases, max_levs);

  // Physics types that won't vary between cases (declare here to avoid calculating the
  // lookup table in the case loop.
  // Sedimentation process.
  Sedimentation sed(constants, true, false);
  // Self-collision processes.
  SelfCollision self_coll;
  // Evaporation process.
  Evaporation evap(constants, sat_form, true, false);

  for (std::size_t icase = 0; icase != num_cases; ++icase) {
    // Grid for this case.
    RainshaftGrid grid = default_grid;
    std::size_t nlev = grid.nlev;

    spaecies::Domain dom;
    spaecies::DimensionPtr lev_dim = dom.add_dimension("level", nlev);
    spaecies::VarDescPtr t_desc = dom.add_var_desc("T", spaecies::Float64Type, {lev_dim}, "K");
    spaecies::VarDescPtr q_desc = dom.add_var_desc("q", spaecies::Float64Type, {lev_dim}, "kg/kg");
    spaecies::VarDescPtr nr_desc = dom.add_var_desc("nr", spaecies::Float64Type, {lev_dim}, "1/kg");
    spaecies::VarDescPtr qr_desc = dom.add_var_desc("qr", spaecies::Float64Type, {lev_dim}, "kg/kg");
    VarDescList state_descs = {t_desc, q_desc, nr_desc, qr_desc};
    VarDescList tend_descs = tend_descs_from_state_descs(dom, state_descs);
    State initial_state(state_descs);

    // Set up initial condition.
    if (initial_condition == "adiabatic") {
      bool converged = warm_adiabatic_initial_condition(constants, grid, sat_form, srf_temp,
                                                        lapse_rate, rel_hum_init, initial_state);
      if (!converged) {
        std::cerr << "Initial condition iteration failed to converge." << std::endl;
        return 1;
      }
    } else {
      throw std::logic_error("Invalid initial condition");
    }
    RainshaftDerivedVars initial_dvars = RainshaftDerivedVars(constants, grid, initial_state);

    // Nudging to initial condition.
    std::shared_ptr<Nudging> nudge;
    // Sum of all processes.
    std::vector<const RainshaftProcess *> exp_process_vec{}, all_process_vec{};
    if (do_nudging) {
      // SPS: Need some kind of span-like interface to avoid having to
      // do this copy.
      VarMut t = initial_state.get_variable("T");
      VarMut q = initial_state.get_variable("q");
      std::vector<double> t_vec, q_vec;
      for (std::size_t i = 0; i != nlev; ++i) {
        t_vec.push_back(t[i]);
        q_vec.push_back(q[i]);
      }
      nudge = std::make_shared<Nudging>(nudge_time_scale, t_vec, q_vec);
      exp_process_vec = {&evap, &*nudge, &self_coll};
      all_process_vec = {&sed, &*nudge, &self_coll, &evap};
    } else {
      exp_process_vec = {&evap, &self_coll};
      all_process_vec = {&sed, &self_coll, &evap};
    }
    SumProcess exp_processes(exp_process_vec);
    // Sum of local processes.
    SumProcess imp_processes = SumProcess{{&sed}};

    SumProcess all_processes = SumProcess{all_process_vec};

    const auto intg = [&]() -> std::unique_ptr<RainshaftIntegrator> {
      if (method_type == "ex") {
        return std::make_unique<ExplicitIntegrator>(constants, grid, &all_processes, state_descs, tend_descs, dt, order, steps_per_output);
      } else if (method_type == "imex") {
        return std::make_unique<IMEXIntegrator>(constants, grid, &exp_processes, &imp_processes, state_descs, tend_descs, dt, order, steps_per_output);
      } else if (method_type == "mri") {
        return std::make_unique<MRIIntegrator>(constants, grid, &imp_processes, &exp_processes, nullptr, state_descs, tend_descs, dt_fast.value(), dt, order, steps_per_output);
      } else {
        throw std::logic_error("Invalid time integration method type");
      }
    }();

    auto before_sol = high_resolution_clock::now();
    RainshaftSolution solution = intg->integrate(initial_time, final_time, initial_state);
    auto after_sol = high_resolution_clock::now();
    // Time taken for solution.
    duration<double, std::milli> walltime_ms = after_sol - before_sol;
    std::cout << "Time: " << walltime_ms.count() << std::endl;
    // Write out grid and all states.
    writer.write_grid(grid, icase);
    writer.write_states(solution.states, icase);
    std::vector<RainshaftDerivedVars> solution_dvars;
    solution_dvars.reserve(solution.states.size());
    for (StateConst state : solution.states) {
      solution_dvars.emplace_back(constants, grid, state);
    }
    writer.write_derived_vars(solution_dvars, icase);
    writer.write_num_rhs_evals(solution.num_rhs_evals, icase);
    writer.write_walltime_ms(walltime_ms.count(), icase);
  }
  return 0;
}
