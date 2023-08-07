#include "spaecies.hpp"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_state.hpp"
#include "rainshaft_tendency.hpp"
#include "rainshaft_solution.hpp"
#include "rainshaft_integrator.hpp"
#include "rainshaft_ncio.hpp"
#include "sedimentation.hpp"
#include "sundials/sundials_context.h"
#include <iostream>

int main(int argc, char** argv)
{
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
  // Time step size in seconds.
  double dt = 1;
  // Time of simulation start.
  double initial_time = 0.;
  // Final time to integrate to.
  double final_time = 300.;
  RainshaftGrid grid = make_e3sm_like_grid(constants, model_top, srf_pres,
                                           srf_temp, lapse_rate);
  sundials::Context sun_ctxt = sundials::Context();
  // Set up initial condition.
  std::vector<double> t, q, nr, qr;
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    t.push_back(srf_temp - lapse_rate * (grid.z_int[il] - 0.5*grid.dz[il]));
    q.push_back(0.);
    nr.push_back(0.);
    qr.push_back(0.);
  }
  RainshaftState initial_state(t, q, nr, qr);
  RainshaftDerivedVars initial_dvars = RainshaftDerivedVars(constants, grid, initial_state);
  // Sedimentation process.
  Sedimentation sed;
  // Evolve state forward.
  RainshaftIntegrator intg(&sun_ctxt, dt);
  RainshaftSolution solution = intg.integrate_ark(sed, initial_time, final_time, constants,
                                                  grid, initial_state, initial_dvars);
  // Write out grid and all states.
  NetcdfWriter writer("./rainshaft_ark2.nc");
  writer.write_grid(grid);
  writer.write_states(solution.states);
  writer.write_derived_vars(solution.dvars);
  writer.write_num_rhs_evals(solution.num_rhs_evals);
  // Ensure that the library is linked and greet the user.
  spaecies::do_nothing();
  return 0;
}
