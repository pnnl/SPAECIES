#include "rainshaft_derived_vars.hpp"
#include <cstddef>
#include <stdexcept>
#include <numeric>

// Calculation of cell widths.
std::vector<double> calc_dz(const RainshaftConstants& constants,
                            const RainshaftGrid& grid,
                            const StateConst& state) {
  // Virtual temperature.
  std::vector<double> t_v(grid.nlev);
  VarConst t = state.get_variable("T").value();
  VarConst q = state.get_variable("q").value();
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    // Virtual temperature factor.
    double t_v_fac = 1. + ((1/constants.epsilon_h2o - 1.) * q[il]);
    t_v[il] = t[il] * t_v_fac;
  }
  return grid.calc_dz(constants, t_v);
}

// Convert cell widths to interface heights.
std::vector<double> dz_to_z_int(const std::vector<double> dz) {
  std::vector<double> z_int(dz.size() + 1);
  z_int.back() = 0.;
  std::partial_sum(dz.crbegin(), dz.crend(), z_int.rbegin() + 1);
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
  VarConst t = state.get_variable("T").value();
  VarConst q = state.get_variable("q").value();
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    rho_dry[il] = rho_dry_from_ideal_gas_law(constants.rdry, constants.epsilon_h2o,
                                             grid.p_mid[il], t[il], q[il]);
  }
  return rho_dry;
}

std::vector<double> calc_lambdar(const RainshaftConstants& constants,
                                 const RainshaftGrid& grid,
                                 const StateConst& state,
                                 const bool regularize) {
  // SPS: If not doing nr limiter here, do it somewhere?
  std::vector<double> lambdar(grid.nlev);
  VarConst nr = state.get_variable("nr").value();
  VarConst qr = state.get_variable("qr").value();
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    double lambdar_value;
    if (regularize) {
      lambdar_value = calc_lambdar(constants.pi, constants.rhow, constants.mur, nr[il], qr[il], constants.qsmall);
    } else {
      lambdar_value = calc_lambdar(constants.pi, constants.rhow, constants.mur, nr[il], qr[il]);
    }
    lambdar[il] = lambdar_value;
  }
  return lambdar;
}

RainshaftDerivedVars::RainshaftDerivedVars(const RainshaftConstants& constants,
                                           const RainshaftGrid& grid,
                                           const StateConst& state,
                                           const bool regularize_lambdar)
  : dz(calc_dz(constants, grid,state)),
    z_int(dz_to_z_int(dz)),
    rho_dry(calc_rho_dry(constants, grid, state)),
    lambdar(calc_lambdar(constants, grid, state, regularize_lambdar)),
    regularize(regularize_lambdar) {
  // Check vector sizes.
  // SPS: Should actually print out mismatched dimensions on failure.
  if (grid.nlev != state.get_variable("T").value().size()) {
    throw std::invalid_argument("grid and state dimensions are mismatched");
  }
}
