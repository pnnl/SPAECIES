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
    const auto breakup = breakup_fac(constants, state.nr[il], state.qr[il]);
    nr_tend[il] = calc_nr_tend(state.nr[il], state.qr[il], breakup, dvars.rho_dry[il]);
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

    const auto breakup = breakup_fac<true>(constants, state.nr[il], state.qr[il]);
    const auto rho_dry = dvars.get_rho_dry<true>(constants, state, il);
    const auto nr_tend = calc_nr_tend<true>(state.nr[il], state.qr[il], breakup, rho_dry);
    const auto [nr_tend_T, nr_tend_q, nr_tend_nr, nr_tend_qr] = get_grad(nr_tend);

    jac(i_nr, i_t) += nr_tend_T;
    jac(i_nr, i_q) += nr_tend_q;
    jac(i_nr, i_nr) += nr_tend_nr;
    jac(i_nr, i_qr) += nr_tend_qr;
  }
}
