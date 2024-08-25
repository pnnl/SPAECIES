#include "self_collision.hpp"
#include "derivatives.hpp"
#include <stdexcept>

RainshaftTendency SelfCollision::calc_tend(const RainshaftConstants &constants,
                                           const RainshaftGrid &grid,
                                           const RainshaftState &state,
                                           const RainshaftDerivedVars &dvars) const
{
  std::vector<double> t_tend(grid.nlev, 0.), q_tend(grid.nlev, 0.);
  std::vector<double> nr_tend(grid.nlev, 0.), qr_tend(grid.nlev, 0.);
  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const auto b = breakup_fac(constants, state.nr[il], state.qr[il]);
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
                                  Matrix jac) const
{
  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const auto i_t = il;
    const auto i_q = i_t + grid.nlev;
    const auto i_nr = i_q + grid.nlev;
    const auto i_qr = i_nr + grid.nlev;

    const auto [b, b_grad] = breakup_fac<true>(constants, state.nr[il], state.qr[il]);
    const auto [db_dnr, db_dqr] = b_grad;

    const auto [dr_dt, dr_dq] = dvars.rho_dry_grad(constants, state, il);

    jac(i_nr, i_t) += -5.78 * b * state.nr[il] * state.qr[il] * dr_dt;
    jac(i_nr, i_q) += -5.78 * b * state.nr[il] * state.qr[il] * dr_dq;
    jac(i_nr, i_nr) += -5.78 * state.qr[il] * dvars.rho_dry[il] * (db_dnr * state.nr[il] + b);
    jac(i_nr, i_qr) += -5.78 * state.nr[il] * dvars.rho_dry[il] * (db_dqr * state.qr[il] + b);
  }
}
