#include "rainshaft_derived_vars.hpp"
#include <cstddef>
#include <cmath>
#include <stdexcept>

// Calculation of cell widths.
std::vector<double> calc_dz(const RainshaftConstants& constants,
                            const RainshaftGrid& grid,
                            const RainshaftState& state) {
  // Virtual temperature.
  std::vector<double> t_v;
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    // Virtual temperature factor.
    double t_v_fac = 1. + ((1/constants.epsilon_h2o - 1.) * state.q[il]);
    t_v.push_back(state.t[il] * t_v_fac);
  }
  return grid.calc_dz(constants, t_v);
}

// Convert cell widths to interface heights.
std::vector<double> dz_to_z_int(const std::vector<double> dz) {
  auto nlev = dz.size();
  std::vector<double> z_int(nlev, 0.);
  // Construct z_int from bottom to top.
  for (int il = nlev-1; il != -1; --il) {
    z_int[il] = z_int[il+1] + dz[il];
  }
  return z_int;
}

// Calculation of dry air density.
// Note that because we are calculating *dry* air density from *total*
// pressure, we use neither the plain temperature nor virtual temperature,
// but instead temperature times (1 + q/epsilon).
std::vector<double> calc_rho_dry(const RainshaftConstants& constants,
                                 const RainshaftGrid& grid,
                                 const RainshaftState& state) {
  std::vector<double> rho_dry;
  for (std::size_t i = 0; i != grid.nlev; ++i) {
    rho_dry.push_back(grid.p_mid[i]
                      / (constants.rdry * state.t[i] *
                         (1 + state.q[i]/constants.epsilon_h2o)));
  }
  return rho_dry;
}

std::vector<double> calc_lambdar(const RainshaftConstants& constants,
                                 const RainshaftGrid& grid,
                                 const RainshaftState& state) {
  // SPS: If not doing nr limiter here, do it somewhere?
  std::vector<double> lambdar;
  for (std::size_t i = 0; i != grid.nlev; ++i) {
    if (state.qr[i] >= constants.qsmall) {
      double lambda_cubed = constants.pi * constants.rhow * state.nr[i]
        / state.qr[i];
      lambdar.push_back(std::cbrt(lambda_cubed));
    } else {
      lambdar.push_back(0.);
    }
  }
  return lambdar;
}

RainshaftDerivedVars::RainshaftDerivedVars(const RainshaftConstants& constants,
                                           const RainshaftGrid& grid,
                                           const RainshaftState& state)
  : dz(calc_dz(constants, grid,state)),
    z_int(dz_to_z_int(dz)),
    rho_dry(calc_rho_dry(constants, grid, state)),
    lambdar(calc_lambdar(constants, grid, state)) {
  // Check vector sizes.
  // SPS: Should actually print out mismatched dimensions on failure.
  if (grid.nlev != state.t.size()) {
    throw std::invalid_argument("grid and state dimensions are mismatched");
  }
}
