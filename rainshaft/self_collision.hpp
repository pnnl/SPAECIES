#ifndef SELF_COLLISION_HPP
#define SELF_COLLISION_HPP
#include <vector>

#include "rainshaft_process.hpp"

class SelfCollision : public RainshaftProcess {

public:

  // Calculate tendency from current state.
  RainshaftTendency calc_tend(const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const RainshaftState& state,
                              const RainshaftDerivedVars& dvars) const;

  RainshaftTendencyJac calc_tend_jac(const RainshaftConstants& constants,
                                           const RainshaftGrid& grid,
                                           const RainshaftState& state,
                                           const RainshaftDerivedVars& dvars) const;

  // For given rain variables, measure whether collisions typically end
  // up merging drops (breakup_fac \approx 1) or is there significant
  // breakup (breakup_fac < 1).
  // Should not be called if nr is near 0.
  double breakup_fac(const RainshaftConstants& constants,
                     double nr, double qr) const;

};

#endif // SELF_COLLISION_HPP
