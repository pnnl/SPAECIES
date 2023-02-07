#ifndef SEDIMENTATION_HPP
#define SEDIMENTATION_HPP
#include <vector>

#include "rainshaft_process.hpp"

class Sedimentation : RainshaftProcess {

public:

  // Calculate tendency from current state.
  RainshaftTendency calc_tend(const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const RainshaftState& state,
                              const RainshaftDerivedVars& dvars);
};

#endif // SEDIMENTATION_HPP
