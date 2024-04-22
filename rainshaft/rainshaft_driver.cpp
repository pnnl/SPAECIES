#include "spaecies.hpp"
#include "evaporation.hpp"
#include "explicit_integrator.hpp"
#include "mri_integrator.hpp"
#include "fixed_substep_integrator.hpp"
#include "forward_euler_integrator.hpp"
#include "nudging.hpp"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_state.hpp"
#include "rainshaft_tendency.hpp"
#include "rainshaft_solution.hpp"
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
#include <math.h>

// // Private function to check function return values
// int check_flag(const int flag, const std::string funcname);

int main(int argc, char** argv)
{
  using std::chrono::high_resolution_clock;
  using std::chrono::duration;

  // Set up model constants.
  // SPS: Choose rho_top in a more principled way?
  RainshaftConstants constants{3.14159265358979323846, // pi
                               287.04, // dry air gas constant
                               1.00464e3, // dry air heat capacity
                               461.50, // water vapor gas constant
                               997., // density of liquid water
                               2.501e6, // latent heat of vaporization
                               0.62197, // Water vapor/dry air molecular mass ratio (unitless)
                               1.e-14, // qsmall
                               1.e-16, // nrsmall
                               9.80616, // gravitational accel.
                               1.e-5, // Minimum allowed mean rain diameter (m)
                               5.e-3, // Maximum allowed mean rain diameter (m)
                               0.988919555598356, // rho at top of column (kg/m^3)
                               1.e3, // nr at top of column (#/kg) (BC)
                               1.e-4}; // qr at top of column (kg/kg) (BC)

  // reusable error flag
  int flag;

  bool uniform_grid_flag = true;

  // Approximate model top in meters.
  // (The grid maker will actually use the next higher-altitude E3SM level.)
  double model_top = 3.e3;
  // Surface pressure and temperature in Pa and K, respectively.
  double srf_pres = 1.e5, srf_temp = 293.15;
  // Lapse rate of initial condition in K/m.
  double lapse_rate = 6.5e-3;
  // Relative humidity of initial condition.
  double rel_hum_init = 0.7;
  // Time scale over which to nudge t and q back to initial condition in seconds.
  double nudge_time_scale = 15. * 60.;
  // Time step size in seconds for substepper i.e. after dt seconds, pass control back to rainshaft code.
  double dt = 1800.00;
  // Time step size in seconds for sundials integrator within each substep
  // double dt_fixedstep = 2.62144;
  double dt_fixedstep = 0.00128 * pow(2, atoi(argv[1]));
  double dt_slowfac = atoi(argv[2]);
  int run = atoi(argv[3]);
  // order of the MRI integrators
  int fastOrder = 3;
  int slowOrder = 3;
  // Time of simulation start.
  double initial_time = 0.;
  // Final time to integrate to.
  double final_time = dt;
  RainshaftGrid grid = make_e3sm_like_grid(constants, model_top, srf_pres,
                                           srf_temp, lapse_rate, uniform_grid_flag);
  auto nlev = grid.nlev;
  // Set up initial condition.
  SaturationFormulae sat_form(constants);
  // Initial condition for nr and qr is just 0, and allocate t and q as well.
  std::vector<double> t(nlev, 0.), q(nlev, 0.), nr(nlev, 0.), qr(nlev, 0.);

  // Coming up with an initial condition for t and q is slightly tricky, because
  // we only have an implicit relationship between t, q, and dz. But since the
  // effect of q on layer height is not large, start by ignoring it, in which
  // case we do have an explicit relationship between t and dz.
  double rog = constants.rdry / constants.g;
  // std::cout << nlev << std::endl;
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
  double t_tolerance = 1.e-1; // Kelvin
  for (int it = 0; it != 1000 && !converged; ++it) {
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

      // std::cout << t[il] << " " << new_temp << std::endl;
      t[il] = new_temp;
    }
  }


  // for (int il = nlev - 1; il != -1; --il) {
  //   std::cout << t[il] << ", " << z_int[il] << ", " << q[il] << std::endl;
  // }

  if (!converged) {
    std::cerr << "Initial condition iteration failed to converge." << std::endl;
    return 1;
  }
  
  RainshaftState initial_state(t, q, nr, qr);
  RainshaftDerivedVars initial_dvars = RainshaftDerivedVars(constants, grid, initial_state);
  // Sedimentation process.
  Sedimentation sed(constants, true, false);

  // Self-collision processes.
  SelfCollision self_coll;
  // Evaporation process.
  Evaporation evap(constants, &sat_form, true, false);
  // Nudging to initial condition.
  Nudging nudge(nudge_time_scale, t, q);

  // Sum of all processes.
  // std::vector<const RainshaftProcess *> micro_processes{&sed, &self_coll, &evap, &nudge};
  std::vector<const RainshaftProcess *> micro_processes{&sed};
  SumProcess all_micro = SumProcess(micro_processes);

  // Sum of local processes.
  std::vector<const RainshaftProcess *> local_processes{&evap, &self_coll, &nudge};
  SumProcess all_local = SumProcess(local_processes);

  std::vector<const RainshaftProcess *> all_processes{&sed, &evap, &self_coll, &nudge};
  SumProcess all_proc = SumProcess(all_processes);

  // Evolve state forward.
  // P3 Settings
  // ForwardEulerIntegrator sed_step(&constants, &grid, &sed, &sun_ctxt);
  // SedCflIntegrator sed_loop(&constants, &grid, &sed, &sed_step);
  // ForwardEulerIntegrator local_step(&constants, &grid, &all_local, &sun_ctxt);
  // std::vector<const RainshaftIntegrator *> seq_ints{&local_step, &sed_loop};
  // SequentialSplitIntegrator seq_step(seq_ints);
  // FixedSubstepIntegrator intg(&seq_step, dt);

  // setup logging
  SUNLogger logger = NULL;

  flag = SUNLogger_Create(
    NULL, // no MPI communicator
    0, // output on process 0 (the only one)
    &logger
  );

  // Attach the logger
  flag = SUNContext_SetLogger(sun_ctxt, logger);

  // Setup log files
  flag = SUNLogger_SetErrorFilename(logger, "stderr");

  flag = SUNLogger_SetWarningFilename(logger, "stderr");

  // logging file name
  flag = SUNLogger_SetInfoFilename(logger, "rainshaft.log");



  // ARKODE Settings
  // ExplicitIntegrator micro_step(&constants, dt_fixedstep, &grid, &all_micro, &sun_ctxt);
  MRIIntegrator micro_step(&constants, fastOrder, slowOrder, dt_fixedstep, dt_slowfac, &grid, &all_local, &all_micro, &all_proc, &sun_ctxt);
  FixedSubstepIntegrator intg(&micro_step, dt);

  // // Pure Forward Euler Settings
  // ForwardEulerIntegrator micro_step(&constants, &grid, &all_micro, &sun_ctxt);
  // FixedSubstepIntegrator intg(&micro_step, dt);
  // Pure Forward Euler Settings
  ForwardEulerIntegrator micro_step(constants, grid, &all_micro);
  FixedSubstepIntegrator intg(&micro_step, dt);
  auto before_sol = high_resolution_clock::now();
  RainshaftSolution solution = intg.integrate(initial_time, final_time, initial_state);
  auto after_sol = high_resolution_clock::now();
  // Time taken for solution.
  duration<double, std::milli> walltime_ms = after_sol - before_sol;
  // Write out grid and all states.

  // file name
  std::string str1 = std::to_string(dt_fixedstep);
  str1.erase(str1.find_first_of('.'), 1);
  str1.erase(0, str1.find_first_not_of('0'));
  str1.erase ( str1.find_last_not_of('0') + 1, std::string::npos );

  std::string str2 = argv[2];
  std::string str3 = argv[3];

  NetcdfWriter writer("./rainshaft_mri_sedfast_t0-1800_order" + std::to_string(fastOrder) + std::to_string(slowOrder) + "_lookup_exp_dt" + str1 + "e-5_outerdt" + str2 + "_run" + str3 + ".nc");
  // NetcdfWriter writer("./rainshaft_t0-1800_order5_exp_dt2e-5_reference.nc");
  // NetcdfWriter writer("./rainshaft_t0-1800_order1_lookup_exp_dt" + str1 + "e-5.nc");
  writer.write_grid(grid);
  writer.write_states(solution.states);
  writer.write_derived_vars(solution.dvars);
  writer.write_num_rhs_evals(solution.num_rhs_evals);
  writer.write_walltime_ms(walltime_ms.count());

  // Ensure that the library is linked and greet the user.
  spaecies::do_nothing();

  if (logger) SUNLogger_Destroy(&logger);  // Free logger

  return 0;
}
