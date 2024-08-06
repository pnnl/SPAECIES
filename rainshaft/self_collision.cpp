#include "self_collision.hpp"
#include "derivatives.hpp"
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include "sunmatrix/sunmatrix_dense.h"
using std::min, std::cbrt, std::exp;

RainshaftTendency SelfCollision::calc_tend(const RainshaftConstants &constants,
                                           const RainshaftGrid &grid,
                                           const RainshaftState &state,
                                           const RainshaftDerivedVars &dvars) const
{
  std::vector<double> t_tend(grid.nlev, 0.), q_tend(grid.nlev, 0.);
  std::vector<double> nr_tend(grid.nlev, 0.), qr_tend(grid.nlev, 0.);
  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const auto [b, b_grad] = breakup_fac(constants, state.nr[il], state.qr[il]);
    nr_tend[il] = -5.78 * b * state.nr[il] * state.qr[il] * dvars.rho_dry[il];
  }
  return RainshaftTendency(t_tend, q_tend, nr_tend, qr_tend);
}

void SelfCollision::calc_tend_jac_prod(const RainshaftConstants &constants,
                                       const RainshaftGrid &grid,
                                       const RainshaftState &state,
                                       const RainshaftDerivedVars &dvars,
                                       const double *const vec,
                                       double *const prod) const
{
  throw std::logic_error("Jacobian product not implemented");
}

void SelfCollision::calc_tend_jac(const RainshaftConstants &constants,
                                  const RainshaftGrid &grid,
                                  const RainshaftState &state,
                                  const RainshaftDerivedVars &dvars,
                                  SUNMatrix jac) const
{
  const auto elem = [jac](const auto i, const auto j) -> auto &
  {
    switch (SUNMatGetID(jac))
    {
    case SUNMATRIX_DENSE:
      return SM_ELEMENT_D(jac, i, j);
    default:
      throw std::logic_error("Unsupported matrix type");
    }
  };

  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const auto i_t = il;
    const auto i_q = i_t + grid.nlev;
    const auto i_nr = i_q + grid.nlev;
    const auto i_qr = i_nr + grid.nlev;

    const auto [b, b_grad] = breakup_fac(constants, state.nr[il], state.qr[il]);
    const auto [db_dnr, db_dqr] = b_grad;

    const auto [dr_dt, dr_dq] = dvars.rho_dry_grad(constants, state, il);

    elem(i_nr, i_t) += -5.78 * b * state.nr[il] * state.qr[il] * dr_dt;
    elem(i_nr, i_q) += -5.78 * b * state.nr[il] * state.qr[il] * dr_dq;
    elem(i_nr, i_nr) += -5.78 * state.qr[il] * dvars.rho_dry[il] * (db_dnr * state.nr[il] + b);
    elem(i_nr, i_qr) += -5.78 * state.nr[il] * dvars.rho_dry[il] * (db_dqr * state.qr[il] + b);
  }
}

ValGrad<2> SelfCollision::breakup_fac(const RainshaftConstants &constants,
                                  double nr, double qr) const
{
  const auto mean_mass_diam = std::cbrt(qr / (constants.pi * constants.rhow * nr));
  const auto breakup_exp_term = std::exp(2300. * (2.8e-4 - mean_mass_diam));
  const auto breakup = 2. - breakup_exp_term;
  if (breakup <= 1) {
    return {breakup, {
      -2300. * breakup_exp_term * mean_mass_diam / (3. * nr),
      2300. * breakup_exp_term * mean_mass_diam / (3. * qr)
    }};
  } else {
    return {1., {0, 0}};
  }
}
