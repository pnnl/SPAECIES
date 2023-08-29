#include "rainshaft_derived_vars.hpp"
#include <cstddef>
#include <cmath>
#include <stdexcept>

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
  : rho_dry(calc_rho_dry(constants, grid, state)),
    lambdar(calc_lambdar(constants, grid, state)) {
  // Check vector sizes.
  // SPS: Should actually print out mismatched dimensions on failure.
  if (grid.nlev != state.t.size()) {
    throw std::invalid_argument("grid and state dimensions are mismatched");
  }
}
