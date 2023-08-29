#ifndef NUDGING_HPP
#define NUDGING_HPP
#include <vector>

#include "rainshaft_process.hpp"

class Nudging : public RainshaftProcess {

public:

  Nudging(double time_scale_in, const std::vector<double>& t,
          const std::vector<double>& q);

  // Calculate tendency from current state.
  RainshaftTendency calc_tend(const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const RainshaftState& state,
                              const RainshaftDerivedVars& dvars) const;

  double time_scale;
  const std::vector<double> t0;
  const std::vector<double> q0;

};

#endif // EVAPORATION_HPP
