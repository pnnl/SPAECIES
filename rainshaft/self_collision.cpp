#include "self_collision.hpp"
#include <cstddef>
#include <cmath>
using std::min, std::cbrt, std::exp;

void SelfCollision::calc_tend(const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const spaecies::State<const double>& state,
                              const RainshaftDerivedVars& dvars,
                              const spaecies::Tendency<double>& tend) const {
  auto nr = state.get_variable("nr");
  auto qr = state.get_variable("qr");
  auto nr_tend = tend.get_variable("nr_tend");
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    double b = breakup_fac(constants, nr[il], qr[il]);
    nr_tend[il] = -5.78 * b * nr[il] * qr[il] * dvars.rho_dry[il];
  }
}

double SelfCollision::breakup_fac(const RainshaftConstants& constants,
                                  double nr, double qr) const {
  // SPS: Should put in an assert that nr is sufficiently large?
  // Diameter of a particle of mean mass (m).
  double mean_mass_diam = std::cbrt(qr / (constants.pi * constants.rhow * nr));
  // Returned factor (dimensionless).
  double b = min(1., 2. - std::exp(2300. * (2.8e-4 - mean_mass_diam)));
  return b;
}
