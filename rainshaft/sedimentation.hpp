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

  // For a given value of lambdar, what are the rain number and mass fall speeds?
  std::vector<double> rain_fall_speeds(const RainshaftConstants& constants,
                                       double rho, double lambdar);

};

#endif // SEDIMENTATION_HPP
