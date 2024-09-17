#ifndef RAINSHAFT_PROCESS_HPP
#define RAINSHAFT_PROCESS_HPP
#include <vector>

#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_derived_vars.hpp"
#include "rainshaft_types.hpp"

class RainshaftProcess {

public:

  // Calculate tendency from current state.
  virtual void calc_tend(const RainshaftConstants& constants,
                         const RainshaftGrid& grid,
                         const StateConst& state,
                         const RainshaftDerivedVars& dvars,
                         const Tendency& tend) const = 0;
};

#endif // RAINSHAFT_PROCESS_HPP
