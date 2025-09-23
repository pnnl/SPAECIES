#ifndef RAINSHAFT_DERIVED_VARS_HPP
#define RAINSHAFT_DERIVED_VARS_HPP
#include <vector>
#include <array>

#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_types.hpp"
#include "derivatives.hpp"

class RainshaftDerivedVars {

public:

  // Constructor from state.
  RainshaftDerivedVars(const RainshaftConstants& constants,
                       const RainshaftGrid& grid,
                       const StateConst& state);

  // Cell heights (m)
  const std::vector<double> dz;
  // Altitude at cell interfaces (m)
  const std::vector<double> z_int;
  // Dry air density (kg/m^3)
  const std::vector<double> rho_dry;
  // Rain size parameter (1/m)
  const std::vector<double> lambdar;

  template <bool WithGrad=false>
  RealOptGrad<WithGrad, 2> get_lambdar(const RainshaftConstants& constants, const double nr, const double qr, std::size_t i) const {
    const double lam = lambdar[i];

    if constexpr (WithGrad) {
      return {lam, {
        lam / (3. * (nr + constants.qsmall * (1.e8))),
        -lam / (3. * (qr + constants.qsmall))
      }};
    } else {
      return lam;
    }
  }

  template <bool WithGrad=false>
  RealOptGrad<WithGrad, 2> get_rho_dry(const RainshaftConstants& constants, const double t, const double q, std::size_t i) const {
    const double rho = rho_dry[i];

    if constexpr (WithGrad) {
      return {rho, {
        -rho / t,
        -rho / (constants.epsilon_h2o + q)
      }};
    } else {
      return rho;
    }
  }

  template <bool WithGrad=false>
  RealOptGrad<WithGrad, 2> get_dz(const RainshaftConstants& constants, const double t, const double q, std::size_t i) const {
    const double dz_i = dz[i];

    if constexpr (WithGrad) {
      const double factor = 1.0/constants.epsilon_h2o - 1.0;
      return {dz_i, {
        dz_i / t,
        dz_i * factor / (1.0 + factor * q)
      }};
    } else {
      return dz_i;
    }
  } 
};

// Convert cell widths to interface heights.
std::vector<double> dz_to_z_int(const std::vector<double> dz);

constexpr double rho_dry_from_ideal_gas_law(double rdry, double epsilon_h2o,
                                            double p_mid, double t, double q) {
  return p_mid / (rdry * t * (1 + q/epsilon_h2o));
}

#endif // RAINSHAFT_DERIVED_VARS_HPP
