#ifndef SELF_COLLISION_HPP
#define SELF_COLLISION_HPP
#include <vector>
#include <cmath>

#include "rainshaft_process.hpp"

class SelfCollision : public RainshaftProcess {

private:
  const bool regularize_lambdar;

  // For given rain variables, measure whether collisions typically end
  // up merging drops (breakup_fac \approx 1) or is there significant
  // breakup (breakup_fac < 1).
  // Should not be called if nr is near 0.
  template<bool WithGrad = false>
  RealOptGrad<WithGrad, 2> breakup_fac(const RainshaftConstants& constants,
                     const double nr, const double qr) const
  {
    double mean_mass_diam, mean_mass_diam_dnr, mean_mass_diam_dqr;
    if (regularize_lambdar) {
      mean_mass_diam = std::cbrt((qr + constants.qsmall) / (constants.pi * constants.rhow * (nr + constants.qsmall * (1.e8))));
      mean_mass_diam_dnr = -mean_mass_diam / (3.0 * (nr + constants.qsmall * (1.e8)));
      mean_mass_diam_dqr = mean_mass_diam / (3.0 * (qr + constants.qsmall));
    } else {
      mean_mass_diam = std::cbrt(qr / (constants.pi * constants.rhow * nr));
      mean_mass_diam_dnr = -mean_mass_diam / (3.0 * nr);
      mean_mass_diam_dqr = mean_mass_diam / (3.0 * qr);
    }
    const auto breakup_exp_term = std::exp(2300. * (mean_mass_diam - 2.8e-4));

    if (constants.epsilon_self_coll > 0.0) {
      // C3 regularization
      const double delta = constants.epsilon_self_coll;
      const auto breakup = (breakup_exp_term >= 1.0 + delta) ? (2.0 - breakup_exp_term) :
                           (breakup_exp_term < 1.0 - delta) ? 1.0 :
                            1.0 - (0.5*std::pow(breakup_exp_term-1.0+delta,2)
                                  - (1.0 + std::cos(M_PI*(breakup_exp_term-1.0)/delta))*std::pow(delta/M_PI,2)
                                  )/(2.0*delta);

      if constexpr (WithGrad) {
        const auto breakup_exp_dnr = 2300. * breakup_exp_term * mean_mass_diam_dnr;
        const auto breakup_exp_dqr = 2300. * breakup_exp_term * mean_mass_diam_dqr;
        const auto dbreakup = -((breakup_exp_term - 1.0 + delta) + std::sin(M_PI*(breakup_exp_term-1.0)/delta) * delta/M_PI) / (2.0*delta);
        const auto breakup_dnr = (breakup_exp_term >= 1.0 + delta) ? (-breakup_exp_dnr) :
                                 (breakup_exp_term < 1.0 - delta) ? 0.0 : dbreakup * breakup_exp_dnr;
        const auto breakup_dqr = (breakup_exp_term >= 1.0 + delta) ? (-breakup_exp_dqr) :
                                 (breakup_exp_term < 1.0 - delta) ? 0.0 : dbreakup * breakup_exp_dqr;
        return {breakup, {breakup_dnr, breakup_dqr}};
      } else {
        return breakup;
      }
    } else {
      // unregularized tendency
      const auto breakup = 2. - breakup_exp_term;
      if (breakup <= 1.0) {
        if constexpr (WithGrad) {
          return {breakup, {
            2300. * breakup_exp_term * mean_mass_diam / (3. * (nr + constants.qsmall * (1.e8))),
            -2300. * breakup_exp_term * mean_mass_diam / (3. * (qr + constants.qsmall))
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

  SelfCollision(const bool regularize_lambdar);

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

  std::set<std::string> get_required_vars() const;
  std::set<std::string> get_optional_vars() const;

};

#endif // SELF_COLLISION_HPP
