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
    * tgamma_lower(3.5, lambdar * 1.3443e-4) * pow(lambdar, -2.5);
  // Integral for 134.43 micron < D <= 1511.64 micron.
  v_evap += sqrt(49.62 * d2g_1third)
    * (tgamma(3., lambdar * 1.3443e-4) - tgamma(3., lambdar * 1.51164e-3))
    / (lambdar * lambdar);
  // Integral for 1511.64 micron < D <= 3477.84 micron.
  v_evap += sqrt(17.32 * d2g_1sixth)
    * (tgamma(2.75, lambdar * 1.51164e-3) - tgamma(2.75, lambdar * 3.47784e-3))
    * pow(lambdar, -1.75);
  // Integral for 3477.84 micron < D.
  v_evap += sqrt(9.17) * tgamma(2.5, lambdar * 3.47784e-3) * pow(lambdar, -1.5);
  return v_evap;
}

// Calculate v_evap using incomplete gamma functions.
double Evaporation::calc_v_evap_numerical(const RainshaftConstants& constants, double lambdar) const {
  double dd = 2.; // micron
  std::size_t low_k = 1;
  std::size_t high_k = 10000;
  double amg_fac = 1000. * constants.pi * constants.rhow / 6.;
  double accum = 0.;
  for (std::size_t kk = low_k; kk != high_k + 1; ++kk) {
    double real_kk = kk;
    double dia_micron = (real_kk - 0.5) * dd;
    double dia = dia_micron * 1.e-6;
    double amg = amg_fac * pow(dia, 3);
    double vt = 0;
    if (dia_micron <= 134.43) {
      vt = 4.5795e3 * pow(amg, 2./3.);
    } else if (dia_micron <= 1511.64) {
      vt = 4.962e1 * pow(amg, 1./3.);
    } else if (dia_micron <= 3477.84) {
      vt = 1.732e1 * pow(amg, 1./6.);
    } else {
      vt = 9.17;
    }
    accum += sqrt(vt * dia) * dia * exp(-lambdar * dia);
  }
  accum *= dd * 1.e-6;
  if (accum < 1.e-30) {
    accum = 1.e-30;
  }
  return accum * lambdar;
}

Evaporation::Evaporation(const RainshaftConstants& constants,
                         const SaturationFormulae* sat_form_in,
                         bool use_v_table,
                         bool use_numerical_integration)
  : sat_form(sat_form_in), use_numerical_integration(use_numerical_integration) {
  if (use_v_table) {
    std::vector<double> range_bounds = {5., 195., 8595.};
    std::vector<double> spacings = {1., 1.};
    std::vector<double> d_microns = LookupTable::calc_x_values(range_bounds,
                                                               spacings);
    std::vector<double> v_values(d_microns.size(), 0.);
    for (std::size_t i = 0; i != d_microns.size(); ++i) {
      double lambdar = 1.e6 / d_microns[i];
      if (use_numerical_integration) {
        v_values[i] = calc_v_evap_numerical(constants, lambdar);
      } else {
        v_values[i] = calc_v_evap_gamma(constants, lambdar);
      }
    }
    v_table.emplace(range_bounds, spacings, v_values);
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

void Evaporation::calc_tend_jac_prod(const RainshaftConstants &constants,
                                 const RainshaftGrid &grid,
                                 const RainshaftState &state,
                                 const RainshaftDerivedVars &dvars,
                                 const double *const vec,
                                 double *const prod) const
{
  throw "Not Implemented";
}

double Evaporation::calc_v_evap(const RainshaftConstants& constants, double lambdar) const {
  if (v_table.has_value()) {
    double d_micron = 1.e6 / lambdar;
    return v_table->lookup_value(d_micron);
  } else {
    if (use_numerical_integration) {
      return calc_v_evap_numerical(constants, lambdar);
    } else {
      return calc_v_evap_gamma(constants, lambdar);
    }
  }
}
