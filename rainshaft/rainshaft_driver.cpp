#include "spaecies.hpp"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_state.hpp"
#include "rainshaft_tendency.hpp"
#include "rainshaft_integrator.hpp"
#include "rainshaft_ncio.hpp"
#include "sedimentation.hpp"
#include <iostream>

int main(int argc, char** argv)
{
  // Set up model constants.
  RainshaftConstants constants{3.14159265358979323846,
                               287.04, 1.00464e3, 461.50, 1000., 2.501e6,
                               0.62197, 1.e-14, 9.80616, 1.e-5, 5.e-3};
  // Approximate model top in meters.
  // (The grid maker will actually use the next higher-altitude E3SM level.)
  double model_top = 2.e3;
  // Surface pressure and temperature in Pa and K, respectively.
  double srf_pres = 1.e5, srf_temp = 293.15;
  // Lapse rate of initial condition in K/m.
  double lapse_rate = 6.5e-3;
  RainshaftGrid grid = make_e3sm_like_grid(constants, model_top, srf_pres,
                                           srf_temp, lapse_rate);
  // Set up initial condition.
  std::vector<double> t, q, nr, qr;
  for (std::size_t i = 0; i != grid.nlev; ++i) {
    t.push_back(srf_temp - lapse_rate * (grid.z_int[i] - 0.5*grid.dz[i]));
    q.push_back(0.);
    nr.push_back(0.);
    qr.push_back(0.);
  }
  RainshaftState initial_state(t, q, nr, qr);
  // Vectors of states and derived variables at different time steps.
  std::vector<RainshaftState> states;
  std::vector<RainshaftDerivedVars> dvars;
  states.push_back(initial_state);
  dvars.push_back(RainshaftDerivedVars(constants, grid, initial_state));
  // Evolve state forward.
  RainshaftIntegrator intg(1.);
  Sedimentation sed;
  RainshaftTendency sed_tend = sed.calc_tend(constants, grid, states[0], dvars[0]);
  RainshaftState new_state = intg.apply_tend(states[0], sed_tend);
  states.push_back(new_state);
  dvars.push_back(RainshaftDerivedVars(constants, grid, new_state));
  // Write out grid and all states.
  NetcdfWriter writer("./rainshaft.nc");
  writer.write_grid(grid);
  writer.write_states(states);
  writer.write_derived_vars(dvars);
  // Ensure that the library is linked and greet the user.
  spaecies::do_nothing();
  return 0;
}
