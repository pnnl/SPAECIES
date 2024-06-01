#include "spaecies.hpp"
#include "evaporation.hpp"
#include "explicit_integrator.hpp"
#include "fixed_substep_integrator.hpp"
// #include "forward_euler_integrator.hpp"
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
#include "imex_integrator.hpp"

int main(int argc, char** argv)
{
  using std::chrono::high_resolution_clock;
  using std::chrono::duration;

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
  // Time step size in seconds.
  double dt = 14;
  // Time of simulation start.
  double initial_time = 0.;
  // Final time to integrate to.
  double final_time = 1800.;
  RainshaftGrid grid = make_e3sm_like_grid(constants, model_top, srf_pres,
                                           srf_temp, lapse_rate);
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
  RainshaftState initial_state(t, q, nr, qr);
  RainshaftDerivedVars initial_dvars = RainshaftDerivedVars(constants, grid, initial_state);
  // Sedimentation process.
  Sedimentation sed(constants, false, false);
  // Self-collision processes.
  SelfCollision self_coll;
  // Evaporation process.
  Evaporation evap(constants, &sat_form, false, false);
  // Nudging to initial condition.
  Nudging nudge(nudge_time_scale, t, q);
  // Sum of all processes.
  SumProcess exp_processes = SumProcess{{&self_coll, &evap, &nudge}};
  // Sum of local processes.
  SumProcess imp_processes = SumProcess{{&sed}};

  SumProcess all_processes = SumProcess{{&sed, &self_coll, &evap, &nudge}};
  // Evolve state forward.
  // P3 Settings
  // ForwardEulerIntegrator sed_step(&constants, &grid, &sed, &sun_ctxt);
  // SedCflIntegrator sed_loop(&constants, &grid, &sed, &sed_step);
  // ForwardEulerIntegrator local_step(&constants, &grid, &all_local, &sun_ctxt);
  // std::vector<const RainshaftIntegrator *> seq_ints{&local_step, &sed_loop};
  // SequentialSplitIntegrator seq_step(seq_ints);
  // FixedSubstepIntegrator intg(&seq_step, dt);
  // ARKODE Settings
  // ExplicitIntegrator micro_step(&constants, &grid, &all_micro, &sun_ctxt);
  // FixedSubstepIntegrator intg(&micro_step, dt);
  // Pure Forward Euler Settings
  // ExplicitIntegrator intg(constants, grid, &all_processes, dt, 2);
  IMEXIntegrator intg(constants, grid, &exp_processes, &imp_processes, dt, 2);
  auto before_sol = high_resolution_clock::now();
  RainshaftSolution solution = intg.integrate(initial_time, final_time, initial_state);
  auto after_sol = high_resolution_clock::now();
  // Time taken for solution.
  duration<double, std::milli> walltime_ms = after_sol - before_sol;
  std::cout << "Time: " << walltime_ms.count() << std::endl;
  // Write out grid and all states.
  NetcdfWriter writer("./rainshaft_1ms_no_table.nc");
  writer.write_grid(grid);
  writer.write_states(solution.states);
  writer.write_derived_vars(solution.dvars);
  writer.write_num_rhs_evals(solution.num_rhs_evals);
  writer.write_walltime_ms(walltime_ms.count());
  // Ensure that the library is linked and greet the user.
  spaecies::do_nothing();
  return 0;
}
