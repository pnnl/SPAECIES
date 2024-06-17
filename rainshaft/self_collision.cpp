#include "self_collision.hpp"
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

    const auto b = breakup_fac(constants, state.nr[il], state.qr[il]);
    elem(i_qr, i_t) += 5.78 * b * state.nr[il] * state.qr[il] * grid.p_mid[il] /
                      (pow(state.t[il], 2) * constants.rdry * (1.0 + state.q[il] / constants.epsilon_h2o));

    elem(i_q, i_qr) += 5.78 * b * state.nr[il] * state.qr[il] * grid.p_mid[il] /
                      (constants.rdry * state.t[il] * constants.epsilon_h2o * pow(1.0 + state.q[il] / constants.epsilon_h2o, 2));

    const auto mean_mass_diam = std::cbrt(state.qr[il] / (constants.pi * constants.rhow * state.nr[il]));
    const auto F = 2. - std::exp(2300. * (2.8e-4 - mean_mass_diam));
    double db;
    if (1.0 <= F)
    {
      db = 0.0;
    }
    else if (dvars.lambdar[il] == 0.0)
    {
      db = -std::exp(2300. * (2.8e-4 - 1.0 / (1.e-5))) * -2300. / (3.0 * (1.e-5));
    }
    else
    {
      db = -std::exp(2300. * (2.8e-4 - 1.0 / dvars.lambdar[il])) * -2300. / (3.0 * dvars.lambdar[il] * state.nr[il]);
    }

    elem(i_nr, i_qr) += -5.78 * state.qr[il] * dvars.rho_dry[il] * (db * state.nr[il] + b);

    elem(i_qr, i_qr) += -5.78 * state.nr[il] * dvars.rho_dry[il] * (-db * state.qr[il] + b);
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
