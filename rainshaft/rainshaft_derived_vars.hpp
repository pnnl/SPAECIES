#ifndef RAINSHAFT_DERIVED_VARS_HPP
#define RAINSHAFT_DERIVED_VARS_HPP
#include <vector>

#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_state.hpp"

class RainshaftDerivedVars {

public:

  // Constructor from state.
  RainshaftDerivedVars(const RainshaftConstants& constants,
                       const RainshaftGrid& grid,
                       const RainshaftState& state);

  // Dry air density (kg/m^3)
  const std::vector<double> rho;
  // Rain size parameter (1/m)
  const std::vector<double> lambdar;
};

#endif // RAINSHAFT_DERIVED_VARS_HPP
