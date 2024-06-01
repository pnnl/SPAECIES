#include "self_collision.hpp"
#include <cstddef>
#include <cmath>
#include <stdexcept>
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
    double b = breakup_fac(constants, state.nr[il], state.qr[il]);
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
  switch (SUNMatGetID(jac))
  {
  default:
    throw std::logic_error("Unsupported matrix type");
  }
}

double SelfCollision::breakup_fac(const RainshaftConstants &constants,
                                  double nr, double qr) const
{
  // SPS: Should put in an assert that nr is sufficiently large?
  // Diameter of a particle of mean mass (m).
  double mean_mass_diam = std::cbrt(qr / (constants.pi * constants.rhow * nr));
  // Returned factor (dimensionless).
  double b = min(1., 2. - std::exp(2300. * (2.8e-4 - mean_mass_diam)));
  return b;
}
