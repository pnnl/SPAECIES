#ifndef RAINSHAFT_DERIVED_VARS_HPP
#define RAINSHAFT_DERIVED_VARS_HPP
#include <vector>

#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_types.hpp"

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
};

// Convert cell widths to interface heights.
std::vector<double> dz_to_z_int(const std::vector<double> dz);

#endif // RAINSHAFT_DERIVED_VARS_HPP
