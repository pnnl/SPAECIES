#include "sedimentation.hpp"
#include <cstddef>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
using std::pow, std::sqrt, std::cbrt, std::exp;
using boost::math::tgamma, boost::math::tgamma_lower;

RainshaftTendency Sedimentation::calc_tend(const RainshaftConstants& constants,
                                           const RainshaftGrid& grid,
                                           const RainshaftState& state,
                                           const RainshaftDerivedVars& dvars) {
  std::vector<double> t_tend, q_tend, nr_tend, qr_tend;
  double nr_tend_lev, qr_tend_lev;
  std::vector<double> vnr, vqr;
  t_tend.push_back(0.);
  q_tend.push_back(0.);
  // SPS: Should have a utility rather than duplicating lambdar calculation.
  double lambdar_top = constants.pi * constants.rhow * constants.nr_top / constants.qr_top;
  lambdar_top = cbrt(lambdar_top);
  std::vector<double> speeds_top = rain_fall_speeds(constants, constants.rho_top, lambdar_top);
  double vnr_top = speeds_top[0], vqr_top = speeds_top[1];
  std::vector<double> speeds = rain_fall_speeds(constants, dvars.rho[0], dvars.lambdar[0]);
  vnr.push_back(speeds[0]);
  vqr.push_back(speeds[1]);
  nr_tend_lev = vnr_top * constants.nr_top * constants.rho_top - vnr[0]*state.nr[0]*dvars.rho[0];
  nr_tend_lev /= grid.dz[0] * dvars.rho[0];
  qr_tend_lev = vqr_top * constants.qr_top * constants.rho_top - vqr[0]*state.qr[0]*dvars.rho[0];
  qr_tend_lev /= grid.dz[0] * dvars.rho[0];
  nr_tend.push_back(nr_tend_lev);
  qr_tend.push_back(qr_tend_lev);
  for (std::size_t il = 1; il != grid.nlev; ++il) {
    t_tend.push_back(0.);
    q_tend.push_back(0.);
    speeds = rain_fall_speeds(constants, dvars.rho[il], dvars.lambdar[il]);
    vnr.push_back(speeds[0]);
    vqr.push_back(speeds[1]);
    nr_tend_lev = vnr[il-1]*state.nr[il-1]*dvars.rho[il-1] - vnr[il]*state.nr[il]*dvars.rho[il];
    nr_tend_lev /= grid.dz[il] * dvars.rho[il];
    qr_tend_lev = vqr[il-1]*state.qr[il-1]*dvars.rho[il-1] - vqr[il]*state.qr[il]*dvars.rho[il];
    qr_tend_lev /= grid.dz[il] * dvars.rho[il];
    nr_tend.push_back(nr_tend_lev);
    qr_tend.push_back(qr_tend_lev);
  }
  return RainshaftTendency(t_tend, q_tend, nr_tend, qr_tend);
}

std::vector<double> Sedimentation::rain_fall_speeds(const RainshaftConstants& constants,
                                                    double rho, double lambdar) {
  // Catches case where there is no rain present.
  if (lambdar == 0.) {
    std::vector<double> speeds = {0., 0.};
    return speeds;
  }
  // Number and mass fall speeds in m/s.
  double vnr(0.), vqr(0.);
  // Factor converting D^3 to drop mass in grams.
  double d3_to_gram = 1000. * constants.pi * constants.rhow / 6.;
  // Gram conversion factor to various powers.
  double d2g_2third = pow(d3_to_gram, 2./3.);
  double d2g_1third = pow(d3_to_gram, 1./3.);
  double d2g_1sixth = pow(d3_to_gram, 1./6.);
  // Start with number term for D <= 134.43 micron.
  vnr += 4579.5 * d2g_2third * tgamma_lower(3, lambdar * 1.3443e-4) / (lambdar * lambdar);
  // Number term for 134.43 micron < D <= 1511.64 micron.
  vnr += 49.62 * d2g_1third
    * (tgamma(2, lambdar * 1.3443e-4) - tgamma(2, lambdar * 1.51164e-3)) / lambdar;
  // Number term for 1511.64 micron < D <= 3477.84 micron.
  vnr += 17.32 * d2g_1sixth
    * (tgamma(1.5, lambdar * 1.51164e-3) - tgamma(1.5, lambdar * 3.47784e-3)) / sqrt(lambdar);
  // Number term for 3477.84 micron < D.
  vnr += 9.17 * exp(-lambdar * 3.47784e-3);
  // Mass term for D <= 134.43 micron.
  vqr += 4579.5 * d2g_2third * tgamma_lower(6, lambdar * 1.3443e-4) / (lambdar * lambdar);
  // Mass term for 134.43 micron < D <= 1511.64 micron.
  vqr += 49.62 * d2g_1third
    * (tgamma(5, lambdar * 1.3443e-4) - tgamma(5, lambdar * 1.51164e-3)) / lambdar;
  // Mass term for 1511.64 micron < D <= 3477.84 micron.
  vqr += 17.32 * d2g_1sixth
    * (tgamma(4.5, lambdar * 1.51164e-3) - tgamma(4.5, lambdar * 3.47784e-3)) / sqrt(lambdar);
  // Mass term for 3477.84 micron < D.
  vqr += 9.17 * tgamma(4, lambdar * 3.47784e-3);
  // Include the normalization gamma(3) for vqr.
  vqr /= 6.;
  // Calculate correction for air density, with reference state at 1000 hPa and 273.15 K.
  double rho_fac = pow(1.e5 / (rho * constants.rdry * 273.15), 0.54);
  std::vector<double> speeds = {rho_fac * vnr, rho_fac * vqr};
  return speeds;
}
