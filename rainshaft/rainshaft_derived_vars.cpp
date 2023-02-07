#include "rainshaft_derived_vars.hpp"
#include <cstddef>
#include <cmath>
#include <stdexcept>

std::vector<double> calc_rho(const RainshaftConstants& constants,
                             const RainshaftGrid& grid,
                             const RainshaftState& state) {
  std::vector<double> rho;
  for (std::size_t i = 0; i != grid.nlev; ++i) {
    rho.push_back(grid.p_mid[i] / (constants.rdry * state.t[i]));
  }
  return rho;
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
  : rho(calc_rho(constants, grid, state)),
    lambdar(calc_lambdar(constants, grid, state)) {
  // Check vector sizes.
  // SPS: Should actually print out mismatched dimensions on failure.
  if (grid.nlev != state.t.size()) {
    throw std::invalid_argument("grid and state dimensions are mismatched");
  }
}
