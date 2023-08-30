#include "evaporation.hpp"
#include "saturation.hpp"
#include <cstddef>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
using std::sqrt, std::pow;
using boost::math::tgamma, boost::math::tgamma_lower;

double calc_diffusivity(double t, double p) {
  return 8.794e-5 * pow(t, 1.81) / p;
}

double calc_dynamic_viscosity(double t) {
  return 1.496e-6 * pow(t, 1.5) / (t + 120.);
}

// Calculate v_evap using incomplete gamma functions.
double Evaporation::calc_v_evap_gamma(const RainshaftConstants& constants, double lambdar) const {
  // Skip function when no rain present.
  if (lambdar == 0.) {
    return 0.;
  }
  double v_evap(0.);
  // Factor converting D^3 to drop mass in grams.
  double d3_to_gram = 1000. * constants.pi * constants.rhow / 6.;
  // Gram conversion factor to various powers.
  double d2g_2third = pow(d3_to_gram, 2./3.);
  double d2g_1third = pow(d3_to_gram, 1./3.);
  double d2g_1sixth = pow(d3_to_gram, 1./6.);
  // Integral for D <= 134.43 micron.
  v_evap = sqrt(4579.5 * d2g_2third)
    * tgamma_lower(17./6., lambdar * 1.3443e-4) * pow(lambdar, -11./6.);
  // Integral for 134.43 micron < D <= 1511.64 micron.
  v_evap += sqrt(49.62 * d2g_1third)
    * (tgamma(8./3., lambdar * 1.3443e-4) - tgamma(8./3., lambdar * 1.51164e-3))
    * pow(lambdar, -5./3.);
  // Integral for 1511.64 micron < D <= 3477.84 micron.
  v_evap += sqrt(17.32 * d2g_1sixth)
    * (tgamma(31./12., lambdar * 1.51164e-3) - tgamma(31./12., lambdar * 3.47784e-3))
    * pow(lambdar, -19./12.);
  // Integral for 3477.84 micron < D.
  v_evap += sqrt(9.17) * tgamma(2.5, lambdar * 3.47784e-3) * pow(lambdar, -1.5);
  return v_evap;
}

Evaporation::Evaporation(const RainshaftConstants& constants,
                         const SaturationFormulae* sat_form_in,
                         bool use_v_table)
  : sat_form(sat_form_in) {
  if (use_v_table) {
    v_table.resize(300);
    for (std::size_t i = 0; i != 20; ++i) {
      double d_micron = 5. + (10. * i);
      double lambdar = 1.e6 / d_micron;
      v_table[i] = calc_v_evap_gamma(constants, lambdar);
    }
    for (std::size_t i = 0; i != 280; ++i) {
      double d_micron = 225. + (30. * i);
      double lambdar = 1.e6 / d_micron;
      v_table[20 + i] = calc_v_evap_gamma(constants, lambdar);
    }
  }
}

RainshaftTendency Evaporation::calc_tend(const RainshaftConstants& constants,
                                         const RainshaftGrid& grid,
                                         const RainshaftState& state,
                                         const RainshaftDerivedVars& dvars) const {
  std::vector<double> t_tend(grid.nlev, 0.), q_tend(grid.nlev, 0.);
  std::vector<double> nr_tend(grid.nlev, 0.), qr_tend(grid.nlev, 0.);
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    // Skip the rest of this if no rain.
    if (state.qr[il] < constants.qsmall) {
      continue;
    }
    double q_sat_dry = sat_form->q_sat_dry(state.t[il], grid.p_mid[il]);
    // Skip the rest of this if not saturated.
    if (q_sat_dry < state.q[il]) {
      continue;
    }
    double dv = calc_diffusivity(state.t[il], grid.p_mid[il]);
    double visc_over_rho = calc_dynamic_viscosity(state.t[il]) / dvars.rho_dry[il];
    double schmidt_num = visc_over_rho / dv;
    double v_evap = calc_v_evap(constants, dvars.lambdar[il]);
    double inv_tau = 2. * constants.pi * state.nr[il] * dvars.rho_dry[il] * dv;
    inv_tau *= 0.78 / dvars.lambdar[il] + (0.32 * pow(schmidt_num, 1./3.)
                                       * sqrt(1./visc_over_rho) * v_evap);
    double abl = 1. + (constants.latvap*constants.latvap * q_sat_dry
                       / (constants.cp * constants.rvapor
                          * state.t[il]*state.t[il]));
    q_tend[il] = (q_sat_dry - state.q[il]) * inv_tau / abl;
    t_tend[il] = - q_tend[il] * (constants.latvap / constants.cp);
    qr_tend[il] = - q_tend[il];
    nr_tend[il] = - q_tend[il] * (state.nr[il] / state.qr[il]);
  }
  return RainshaftTendency(t_tend, q_tend, nr_tend, qr_tend);
}

double Evaporation::calc_v_evap(const RainshaftConstants& constants, double lambdar) const {
  if (v_table.size() > 0) {
    double d_micron = 1.e6 / lambdar;
    if (d_micron < 5.) {
      return v_table[0];
    } else if (d_micron < 195.) {
      std::size_t low_ind = (((std::size_t) d_micron) - 5) / 10;
      double frac_part = d_micron - ((low_ind * 10) + 5);
      return (1. - frac_part) * v_table[low_ind]
        + frac_part * v_table[low_ind+1];
    } else if (d_micron < 8595.) {
      std::size_t low_ind = (((std::size_t) d_micron - 195)) / 30 + 19;
      double frac_part = d_micron - (((low_ind - 19) * 30) + 195);
      return (1. - frac_part) * v_table[low_ind]
        + frac_part * v_table[low_ind+1];
    } else {
      return v_table[299];
    }
  } else {
    return calc_v_evap_gamma(constants, lambdar);
  }
}
