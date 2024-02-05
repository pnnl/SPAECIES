#ifndef RAINSHAFT_PROCESS_HPP
#define RAINSHAFT_PROCESS_HPP
#include <vector>

#include "spaecies.hpp"

#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_derived_vars.hpp"
#include "rainshaft_tendency.hpp"

class RainshaftProcess {

public:

  // Calculate tendency from current state.
  virtual RainshaftTendency calc_tend(const RainshaftConstants& constants,
                                      const RainshaftGrid& grid,
                                      const spaecies::VariableArrayView<double>& state,
                                      const RainshaftDerivedVars& dvars) const = 0;
};

#endif // RAINSHAFT_PROCESS_HPP
