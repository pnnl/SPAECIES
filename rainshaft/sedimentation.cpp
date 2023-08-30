#include "sedimentation.hpp"
#include <cstddef>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
using std::pow, std::sqrt, std::cbrt, std::exp;
using boost::math::tgamma, boost::math::tgamma_lower;

Sedimentation::Sedimentation(const RainshaftConstants& constants, bool use_v_table) {
  if (use_v_table) {
    v0_table.resize(300);
    v3_table.resize(300);
    for (std::size_t i = 0; i != 20; ++i) {
      double d_micron = 5. + (10. * i);
      double lambdar = 1.e6 / d_micron;
      auto speeds = rain_fall_speeds_stp_gamma(constants, lambdar);
      v0_table[i] = speeds[0];
      v3_table[i] = speeds[1];
    }
    for (std::size_t i = 0; i != 280; ++i) {
      double d_micron = 225. + (30. * i);
      double lambdar = 1.e6 / d_micron;
      auto speeds = rain_fall_speeds_stp_gamma(constants, lambdar);
      v0_table[20 + i] = speeds[0];
      v3_table[20 + i] = speeds[1];
    }
  }
}

RainshaftTendency Sedimentation::calc_tend(const RainshaftConstants& constants,
                                           const RainshaftGrid& grid,
                                           const RainshaftState& state,
                                           const RainshaftDerivedVars& dvars) const {
  std::vector<double> t_tend, q_tend, nr_tend, qr_tend;
  double nr_tend_lev, qr_tend_lev;
  std::vector<double> v0, v3;
  t_tend.push_back(0.);
  q_tend.push_back(0.);
  // SPS: Should have a utility rather than duplicating lambdar calculation.
  double lambdar_top = constants.pi * constants.rhow * constants.nr_top / constants.qr_top;
  lambdar_top = cbrt(lambdar_top);
  std::vector<double> speeds_top = rain_fall_speeds(constants, constants.rho_top, lambdar_top);
  double v0_top = speeds_top[0], v3_top = speeds_top[1];
  std::vector<double> speeds = rain_fall_speeds(constants, dvars.rho_dry[0], dvars.lambdar[0]);
  v0.push_back(speeds[0]);
  v3.push_back(speeds[1]);
  nr_tend_lev = v0_top * constants.nr_top * constants.rho_top - v0[0]*state.nr[0]*dvars.rho_dry[0];
  nr_tend_lev /= dvars.dz[0] * dvars.rho_dry[0];
  qr_tend_lev = v3_top * constants.qr_top * constants.rho_top - v3[0]*state.qr[0]*dvars.rho_dry[0];
  qr_tend_lev /= dvars.dz[0] * dvars.rho_dry[0];
  nr_tend.push_back(nr_tend_lev);
  qr_tend.push_back(qr_tend_lev);
  for (std::size_t il = 1; il != grid.nlev; ++il) {
    t_tend.push_back(0.);
    q_tend.push_back(0.);
    speeds = rain_fall_speeds(constants, dvars.rho_dry[il], dvars.lambdar[il]);
    v0.push_back(speeds[0]);
    v3.push_back(speeds[1]);
    nr_tend_lev = v0[il-1]*state.nr[il-1]*dvars.rho_dry[il-1] - v0[il]*state.nr[il]*dvars.rho_dry[il];
    nr_tend_lev /= dvars.dz[il] * dvars.rho_dry[il];
    qr_tend_lev = v3[il-1]*state.qr[il-1]*dvars.rho_dry[il-1] - v3[il]*state.qr[il]*dvars.rho_dry[il];
    qr_tend_lev /= dvars.dz[il] * dvars.rho_dry[il];
    nr_tend.push_back(nr_tend_lev);
    qr_tend.push_back(qr_tend_lev);
  }
  return RainshaftTendency(t_tend, q_tend, nr_tend, qr_tend);
}

std::vector<double> Sedimentation::rain_fall_speeds(const RainshaftConstants& constants,
                                                    double rho_dry, double lambdar) const {
  // Catches case where there is no rain present.
  if (lambdar == 0.) {
    std::vector<double> speeds = {0., 0.};
    return speeds;
  }
  auto speeds = rain_fall_speeds_stp(constants, lambdar);
  // Calculate correction for air density, with reference state at 1000 hPa and 273.15 K.
  double rf = rho_fac(constants, rho_dry);
  speeds[0] *= rf;
  speeds[1] *= rf;
  return speeds;
}

std::vector<double> Sedimentation::rain_fall_speeds_stp(const RainshaftConstants& constants,
                                                        double lambdar) const {
  if (v0_table.size() > 0) {
    double d_micron = 1.e6 / lambdar;
    double v0, v3;
    if (d_micron < 5.) {
      v0 = v0_table[0];
      v3 = v3_table[0];
    } else if (d_micron < 195.) {
      std::size_t low_ind = (((std::size_t) d_micron) - 5) / 10;
      double frac_part = d_micron - ((low_ind * 10) + 5);
      v0 = (1. - frac_part) * v0_table[low_ind] + frac_part * v0_table[low_ind+1];
      v3 = (1. - frac_part) * v3_table[low_ind] + frac_part * v3_table[low_ind+1];
    } else if (d_micron < 8595.) {
      std::size_t low_ind = (((std::size_t) d_micron - 195)) / 30 + 19;
      double frac_part = d_micron - (((low_ind - 19) * 30) + 195);
      v0 = (1. - frac_part) * v0_table[low_ind] + frac_part * v0_table[low_ind+1];
      v3 = (1. - frac_part) * v3_table[low_ind] + frac_part * v3_table[low_ind+1];
    } else {
      v0 = v0_table[299];
      v3 = v3_table[299];
    }
    return std::vector<double>{v0, v3};
  } else {
    return rain_fall_speeds_stp_gamma(constants, lambdar);
  }
}

std::vector<double> Sedimentation::rain_fall_speeds_stp_gamma(const RainshaftConstants& constants,
                                                              double lambdar) const {
  // Catches case where there is no rain present.
  if (lambdar == 0.) {
    std::vector<double> speeds = {0., 0.};
    return speeds;
  }
  // Number and mass fall speeds in m/s.
  double v0(0.), v3(0.);
  // Factor converting D^3 to drop mass in grams.
  double d3_to_gram = 1000. * constants.pi * constants.rhow / 6.;
  // Gram conversion factor to various powers.
  double d2g_2third = pow(d3_to_gram, 2./3.);
  double d2g_1third = pow(d3_to_gram, 1./3.);
  double d2g_1sixth = pow(d3_to_gram, 1./6.);
  // Start with number term for D <= 134.43 micron.
  v0 += 4579.5 * d2g_2third * tgamma_lower(3, lambdar * 1.3443e-4) / (lambdar * lambdar);
  // Number term for 134.43 micron < D <= 1511.64 micron.
  v0 += 49.62 * d2g_1third
    * (tgamma(2, lambdar * 1.3443e-4) - tgamma(2, lambdar * 1.51164e-3)) / lambdar;
  // Number term for 1511.64 micron < D <= 3477.84 micron.
  v0 += 17.32 * d2g_1sixth
    * (tgamma(1.5, lambdar * 1.51164e-3) - tgamma(1.5, lambdar * 3.47784e-3)) / sqrt(lambdar);
  // Number term for 3477.84 micron < D.
  v0 += 9.17 * exp(-lambdar * 3.47784e-3);
  // Mass term for D <= 134.43 micron.
  v3 += 4579.5 * d2g_2third * tgamma_lower(6, lambdar * 1.3443e-4) / (lambdar * lambdar);
  // Mass term for 134.43 micron < D <= 1511.64 micron.
  v3 += 49.62 * d2g_1third
    * (tgamma(5, lambdar * 1.3443e-4) - tgamma(5, lambdar * 1.51164e-3)) / lambdar;
  // Mass term for 1511.64 micron < D <= 3477.84 micron.
  v3 += 17.32 * d2g_1sixth
    * (tgamma(4.5, lambdar * 1.51164e-3) - tgamma(4.5, lambdar * 3.47784e-3)) / sqrt(lambdar);
  // Mass term for 3477.84 micron < D.
  v3 += 9.17 * tgamma(4, lambdar * 3.47784e-3);
  // Include the normalization gamma(3) for v3.
  v3 /= 6.;
  std::vector<double> speeds = {v0, v3};
  return speeds;
}

double Sedimentation::rho_fac(const RainshaftConstants& constants, double rho_dry) const {
  return pow(1.e5 / (rho_dry * constants.rdry * 273.15), 0.54);
}
