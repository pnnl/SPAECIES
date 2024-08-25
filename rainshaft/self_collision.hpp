#ifndef SELF_COLLISION_HPP
#define SELF_COLLISION_HPP
#include <vector>
#include <cmath>

#include "rainshaft_process.hpp"

class SelfCollision : public RainshaftProcess {

private:
  // For given rain variables, measure whether collisions typically end
  // up merging drops (breakup_fac \approx 1) or is there significant
  // breakup (breakup_fac < 1).
  // Should not be called if nr is near 0.
  template<bool WithGrad = false>
  auto breakup_fac(const RainshaftConstants& constants,
                     const double nr, const double qr) const
  {
    const auto mean_mass_diam = std::cbrt(qr / (constants.pi * constants.rhow * nr));
    const auto breakup_exp_term = std::exp(2300. * (2.8e-4 - mean_mass_diam));
    const auto breakup = 2. - breakup_exp_term;

    if constexpr (WithGrad) {
      if (breakup <= 1.0) {
        return ValGrad<2>{
          breakup,
          {
            -2300. * breakup_exp_term * mean_mass_diam / (3. * nr),
            2300. * breakup_exp_term * mean_mass_diam / (3. * qr)
          }
        };
      } else {
        return ValGrad<2>{1.0, {0.0, 0.0}};
      }
    } else {
      return std::min(1.0, breakup);
    }
  }


public:

  // Calculate tendency from current state.
  RainshaftTendency calc_tend(const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const RainshaftState& state,
                              const RainshaftDerivedVars& dvars) const;

  void calc_tend_jac_prod(const RainshaftConstants &constants,
                          const RainshaftGrid &grid,
                          const RainshaftState &state,
                          const RainshaftDerivedVars &dvars,
                          const double *const vec,
                          double *const prod) const;

  void calc_tend_jac(const RainshaftConstants &constants,
                             const RainshaftGrid &grid,
                             const RainshaftState &state,
                             const RainshaftDerivedVars &dvars,
                             Matrix jac) const;

};

#endif // SELF_COLLISION_HPP
