#ifndef SEDIMENTATION_HPP
#define SEDIMENTATION_HPP
#include <vector>

#include "rainshaft_process.hpp"

class Sedimentation : public RainshaftProcess {

public:

  // Calculate tendency from current state.
  RainshaftTendency calc_tend(const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const RainshaftState& state,
                              const RainshaftDerivedVars& dvars) const;

  // For a given value of lambdar, what are the rain number and mass fall speeds?
  std::vector<double> rain_fall_speeds(const RainshaftConstants& constants,
                                       double rho, double lambdar) const;

};

#endif // SEDIMENTATION_HPP
