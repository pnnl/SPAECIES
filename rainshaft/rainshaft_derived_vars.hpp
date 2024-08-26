#ifndef RAINSHAFT_DERIVED_VARS_HPP
#define RAINSHAFT_DERIVED_VARS_HPP
#include <vector>
#include <array>

#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_state.hpp"
#include "derivatives.hpp"

class RainshaftDerivedVars {

public:

  // Constructor from state.
  RainshaftDerivedVars(const RainshaftConstants& constants,
                       const RainshaftGrid& grid,
                       const RainshaftState& state);

  // Cell heights (m)
  const std::vector<double> dz;
  // Altitude at cell interfaces (m)
  const std::vector<double> z_int;
  // Dry air density (kg/m^3)
  const std::vector<double> rho_dry;
  // Rain size parameter (1/m)
  const std::vector<double> lambdar;

  template <bool WithGrad=false>
  Val<WithGrad, 2> get_lambdar(const RainshaftConstants& constants, const RainshaftState& state, std::size_t i) const {
    const auto lam = lambdar[i];

    if constexpr (WithGrad) {
      return {lam, {
        state.qr[i] >= constants.qsmall ? lam / (3. * state.nr[i]) : 0.0,
        state.qr[i] >= constants.qsmall ? -lam / (3. * state.qr[i]) : 0.0
      }};
    } else {
      return lam;
    }
  }

  template <bool WithGrad=false>
  Val<WithGrad, 2> get_rho_dry(const RainshaftConstants& constants, const RainshaftState& state, std::size_t i) const {
    const auto rho = rho_dry[i];

    if constexpr (WithGrad) {
      return {rho, {
        -rho / state.t[i],
        -rho / (constants.epsilon_h2o + state.q[i])
      }};
    } else {
      return rho;
    }
  }
};

// Convert cell widths to interface heights.
std::vector<double> dz_to_z_int(const std::vector<double> dz);

#endif // RAINSHAFT_DERIVED_VARS_HPP
