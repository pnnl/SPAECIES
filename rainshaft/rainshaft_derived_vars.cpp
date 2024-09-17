#include "rainshaft_derived_vars.hpp"
#include <cstddef>
#include <cmath>
#include <stdexcept>

// Calculation of cell widths.
std::vector<double> calc_dz(const RainshaftConstants& constants,
                            const RainshaftGrid& grid,
                            const StateConst& state) {
  // Virtual temperature.
  std::vector<double> t_v(grid.nlev);
  VarConst t = state.get_variable("T");
  VarConst q = state.get_variable("q");
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    // Virtual temperature factor.
    double t_v_fac = 1. + ((1/constants.epsilon_h2o - 1.) * q[il]);
    t_v[il] = t[il] * t_v_fac;
  }
  return grid.calc_dz(constants, t_v);
}

// Convert cell widths to interface heights.
std::vector<double> dz_to_z_int(const std::vector<double> dz) {
  std::size_t nlev = dz.size();
  std::vector<double> z_int(nlev+1, 0.);
  // Construct z_int from bottom to top.
  for (std::size_t il = nlev-1; il != -1; --il) {
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
                                 const StateConst& state) {
  std::vector<double> rho_dry(grid.nlev);
  VarConst t = state.get_variable("T");
  VarConst q = state.get_variable("q");
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    rho_dry[il] = grid.p_mid[il] / (constants.rdry * t[il] *
                                  (1 + q[il]/constants.epsilon_h2o));
  }
  return rho_dry;
}

std::vector<double> calc_lambdar(const RainshaftConstants& constants,
                                 const RainshaftGrid& grid,
                                 const StateConst& state) {
  // SPS: If not doing nr limiter here, do it somewhere?
  std::vector<double> lambdar(grid.nlev);
  VarConst nr = state.get_variable("nr");
  VarConst qr = state.get_variable("qr");
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    if (qr[il] >= constants.qsmall) {
      double lambda_cubed = constants.pi * constants.rhow * nr[il]
        / qr[il];
      lambdar[il] = std::cbrt(lambda_cubed);
    }
  }
  return lambdar;
}

RainshaftDerivedVars::RainshaftDerivedVars(const RainshaftConstants& constants,
                                           const RainshaftGrid& grid,
                                           const StateConst& state)
  : dz(calc_dz(constants, grid,state)),
    z_int(dz_to_z_int(dz)),
    rho_dry(calc_rho_dry(constants, grid, state)),
    lambdar(calc_lambdar(constants, grid, state)) {
  // Check vector sizes.
  // SPS: Should actually print out mismatched dimensions on failure.
  if (grid.nlev != state.get_variable("T").size()) {
    throw std::invalid_argument("grid and state dimensions are mismatched");
  }
}
