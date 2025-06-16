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
  RealOptGrad<WithGrad, 2> breakup_fac(const RainshaftConstants& constants,
                     const double nr, const double qr) const
  {
    const auto mean_mass_diam = std::cbrt(qr / (constants.pi * constants.rhow * nr));
    const auto breakup_exp_term = std::exp(2300. * (mean_mass_diam - 2.8e-4));
    const auto breakup = 2. - breakup_exp_term;

    if (breakup <= 1.0) {
      if constexpr (WithGrad) {
        return {breakup, {
          2300. * breakup_exp_term * mean_mass_diam / (3. * nr),
          -2300. * breakup_exp_term * mean_mass_diam / (3. * qr)
        }};
      } else {
        return breakup;
      }
    } else {
      if constexpr (WithGrad) {
        return {1.0, {0.0, 0.0}};
      } else {
        return 1.0;
      }
    }
  }

  template<bool WithGrad = false>
  RealOptGrad<WithGrad, 4> calc_nr_tend(const double nr,
                                const double qr,
                                const RealOptGrad<WithGrad, 2> breakup,
                                const RealOptGrad<WithGrad, 2> rho_dry) const
  {
    const auto nr_tend_fac =  -5.78 * get_val(breakup) * nr * qr;
    const auto nr_tend = nr_tend_fac * get_val(rho_dry);

    if constexpr (WithGrad) {
      const auto [rho_dry_dT, rho_dry_dq] = get_grad(rho_dry);
      const auto [breakup_dnr, breakup_dqr] = get_grad(breakup);
      return {nr_tend, {
        nr_tend_fac * rho_dry_dT,
        nr_tend_fac * rho_dry_dq,
        -5.78 * qr * get_val(rho_dry) * (breakup_dnr * nr + get_val(breakup)),
        -5.78 * nr * get_val(rho_dry) * (breakup_dqr * qr + get_val(breakup))
      }};
    } else {
      return nr_tend;
    }
  }


public:

  // Calculate tendency from current state.
  void calc_tend(const RainshaftConstants& constants,
                 const RainshaftGrid& grid,
                 const StateConst& state,
                 const RainshaftDerivedVars& dvars,
                 Tendency& tend) const;

  void calc_tend_jac(const RainshaftConstants &constants,
                             const RainshaftGrid &grid,
                             const StateConst& state,
                             const RainshaftDerivedVars &dvars,
                             Matrix jac) const;

};

#endif // SELF_COLLISION_HPP
