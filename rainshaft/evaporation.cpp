#include "evaporation.hpp"
#include "saturation.hpp"
#include <cstddef>
#include <cmath>
#include "sunmatrix/sunmatrix_dense.h"
#include <boost/math/special_functions/gamma.hpp>
using boost::math::tgamma, boost::math::tgamma_lower;
using std::sqrt, std::pow;

double calc_diffusivity(double t, double p) {
  return 8.794e-5 * pow(t, 1.81) / p;
}

double calc_diffusivity_dT(double t, double p) {
  return 8.794e-5 * 1.81 / p * pow(t, 0.81);
}

double calc_dynamic_viscosity(double t) {
  return 1.496e-6 * pow(t, 1.5) / (t + 120.);
}

double calc_dynamic_viscosity_dT(double t) {
  return 1.496e-6 * ((t + 120.)*1.5*pow(t, 0.5) - pow(t, 1.5)) / pow(t + 120., 2);
}

// Calculate v_evap using incomplete gamma functions.
double Evaporation::calc_v_evap_gamma(const RainshaftConstants &constants, double lambdar) const
{
  // Skip function when no rain present.
  if (lambdar == 0.)
  {
    return 0.;
  }
  double v_evap(0.);
  // Factor converting D^3 to drop mass in grams.
  double d3_to_gram = 1000. * constants.pi * constants.rhow / 6.;
  // Gram conversion factor to various powers.
  double d2g_2third = pow(d3_to_gram, 2. / 3.);
  double d2g_1third = pow(d3_to_gram, 1. / 3.);
  double d2g_1sixth = pow(d3_to_gram, 1. / 6.);
  // Integral for D <= 134.43 micron.
  v_evap = sqrt(4579.5 * d2g_2third) * tgamma_lower(3.5, lambdar * 1.3443e-4) * pow(lambdar, -2.5);
  // Integral for 134.43 micron < D <= 1511.64 micron.
  v_evap += sqrt(49.62 * d2g_1third) * (tgamma(3., lambdar * 1.3443e-4) - tgamma(3., lambdar * 1.51164e-3)) / (lambdar * lambdar);
  // Integral for 1511.64 micron < D <= 3477.84 micron.
  v_evap += sqrt(17.32 * d2g_1sixth) * (tgamma(2.75, lambdar * 1.51164e-3) - tgamma(2.75, lambdar * 3.47784e-3)) * pow(lambdar, -1.75);
  // Integral for 3477.84 micron < D.
  v_evap += sqrt(9.17) * tgamma(2.5, lambdar * 3.47784e-3) * pow(lambdar, -1.5);
  return v_evap;
}

// Calculate v_evap using incomplete gamma functions.
double Evaporation::calc_v_evap_numerical(const RainshaftConstants &constants, double lambdar) const
{
  double dd = 2.; // micron
  std::size_t low_k = 1;
  std::size_t high_k = 10000;
  double amg_fac = 1000. * constants.pi * constants.rhow / 6.;
  double accum = 0.;
  for (std::size_t kk = low_k; kk != high_k + 1; ++kk)
  {
    double real_kk = kk;
    double dia_micron = (real_kk - 0.5) * dd;
    double dia = dia_micron * 1.e-6;
    double amg = amg_fac * pow(dia, 3);
    double vt = 0;
    if (dia_micron <= 134.43)
    {
      vt = 4.5795e3 * pow(amg, 2. / 3.);
    }
    else if (dia_micron <= 1511.64)
    {
      vt = 4.962e1 * pow(amg, 1. / 3.);
    }
    else if (dia_micron <= 3477.84)
    {
      vt = 1.732e1 * pow(amg, 1. / 6.);
    }
    else
    {
      vt = 9.17;
    }
    accum += sqrt(vt * dia) * dia * exp(-lambdar * dia);
  }
  accum *= dd * 1.e-6;
  if (accum < 1.e-30)
  {
    accum = 1.e-30;
  }
  return accum * lambdar;
}

Evaporation::Evaporation(const RainshaftConstants &constants,
                         const SaturationFormulae *sat_form_in,
                         bool use_v_table,
                         bool use_numerical_integration)
    : sat_form(sat_form_in), use_numerical_integration(use_numerical_integration)
{
  if (use_v_table)
  {
    std::vector<double> range_bounds = {5., 195., 8595.};
    std::vector<double> spacings = {1., 1.};
    std::vector<double> d_microns = LookupTable::calc_x_values(range_bounds,
                                                               spacings);
    std::vector<double> v_values(d_microns.size(), 0.);
    std::vector<double> dv_values(d_microns.size(), 0.);
    for (std::size_t i = 0; i != d_microns.size(); ++i)
    {
      double lambdar = 1.e6 / d_microns[i];
      if (use_numerical_integration)
      {
        v_values[i] = calc_v_evap_numerical(constants, lambdar);
      }
      else
      {
        v_values[i] = calc_v_evap_gamma(constants, lambdar);
        dv_values[i] = calc_v_evap_dlambda_gamma(constants, lambdar);
      }
    }
    v_table.emplace(range_bounds, spacings, v_values);
    dv_table.emplace(range_bounds, spacings, dv_values);
  }
}

RainshaftTendency Evaporation::calc_tend(const RainshaftConstants &constants,
                                         const RainshaftGrid &grid,
                                         const RainshaftState &state,
                                         const RainshaftDerivedVars &dvars) const
{
  std::vector<double> t_tend(grid.nlev, 0.), q_tend(grid.nlev, 0.);
  std::vector<double> nr_tend(grid.nlev, 0.), qr_tend(grid.nlev, 0.);
  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    // Skip the rest of this if no rain.
    if (state.qr[il] < constants.qsmall)
    {
      continue;
    }
    double q_sat_dry = sat_form->q_sat_dry(state.t[il], grid.p_mid[il]);
    // Skip the rest of this if not saturated.
    if (q_sat_dry < state.q[il])
    {
      continue;
    }
    double dv = calc_diffusivity(state.t[il], grid.p_mid[il]);
    double visc_over_rho = calc_dynamic_viscosity(state.t[il]) / dvars.rho_dry[il];
    double schmidt_num = visc_over_rho / dv;
    double v_evap = calc_v_evap(constants, dvars.lambdar[il]);
    double inv_tau = 2. * constants.pi * state.nr[il] * dvars.rho_dry[il] * dv;
    inv_tau *= 0.78 / dvars.lambdar[il] + (0.32 * pow(schmidt_num, 1. / 3.) * sqrt(1. / visc_over_rho) * v_evap);
    double abl = 1. + (constants.latvap * constants.latvap * q_sat_dry / (constants.cp * constants.rvapor * state.t[il] * state.t[il]));
    q_tend[il] = (q_sat_dry - state.q[il]) * inv_tau / abl;
    t_tend[il] = -q_tend[il] * (constants.latvap / constants.cp);
    qr_tend[il] = -q_tend[il];
    nr_tend[il] = -q_tend[il] * (state.nr[il] / state.qr[il]);
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
  throw std::logic_error("Jacobian product not implemented");
}

void Evaporation::calc_tend_jac(const RainshaftConstants &constants,
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
    const auto q_sat_dry = sat_form->q_sat_dry(state.t[il], grid.p_mid[il]);

    if (state.qr[il] < constants.qsmall || q_sat_dry < state.q[il])
    {
      continue;
    }

    const auto q_sat_dry_dT = sat_form->q_sat_dry_dT(state.t[il], grid.p_mid[il]);

    const auto dv = calc_diffusivity(state.t[il], grid.p_mid[il]);
    const auto dv_dT = calc_diffusivity_dT(state.t[il], grid.p_mid[il]);

    const auto rho_dry_dT = dvars.rho_dry[il] * -1.0 / state.t[il];
    const auto rho_dry_dq = dvars.rho_dry[il] * -1.0 / (constants.epsilon_h2o + state.q[il]);
    const auto visc_over_rho = calc_dynamic_viscosity(state.t[il]) / dvars.rho_dry[il];
    const auto visc_over_rho_dT = (dvars.rho_dry[il] * calc_dynamic_viscosity_dT(state.t[il]) - calc_dynamic_viscosity(state.t[il]) * rho_dry_dT) / pow(dvars.rho_dry[il], 2);
    const auto visc_over_rho_dq = -calc_dynamic_viscosity(state.t[il]) / pow(dvars.rho_dry[il], 2) * rho_dry_dq;

    const auto schmidt_num = visc_over_rho / dv;
    const auto schmidt_num_dT = (dv * visc_over_rho_dT - visc_over_rho * dv_dT) / pow(dv, 2);

    // evaporation fall speed
    const auto v_evap = calc_v_evap(constants, dvars.lambdar[il]);
    const auto v_evap_dlambda = calc_v_evap_dlambda(constants, dvars.lambdar[il]);

    // psychometric correction factor abl
    const auto abl = 1. + (constants.latvap * constants.latvap * q_sat_dry / (constants.cp * constants.rvapor * state.t[il] * state.t[il]));
    const auto abl_dT = constants.latvap * constants.latvap / (constants.cp * constants.rvapor) *
                        (pow(state.t[il], 2) * q_sat_dry_dT - 2.0 * q_sat_dry * state.t[il]) / pow(state.t[il], 4);

    const auto lambdar_dnr = dvars.lambdar[il] / (3.0 * state.nr[il]);
    const auto lambdar_dqr = -dvars.lambdar[il] / (3.0 * state.qr[il]);

    const auto tau_fac = 0.78 / dvars.lambdar[il] + (0.32 * pow(schmidt_num, 1. / 3.) * sqrt(1. / visc_over_rho) * v_evap);
    const auto tau_fac_dT = 0.32 * v_evap * ((1.0 / 3.0) * pow(schmidt_num, -2. / 3.) * schmidt_num_dT * sqrt(1. / visc_over_rho) + pow(schmidt_num, 1. / 3.) * 0.5 / sqrt(1. / visc_over_rho) * -1.0 / pow(visc_over_rho, 2) * visc_over_rho_dT);
    const auto tau_fac_dq = 0.32 * pow(schmidt_num, 1. / 3.) * v_evap * 0.5 / sqrt(1. / visc_over_rho) *
                            -1.0 / pow(visc_over_rho, 2) * visc_over_rho_dq;
    const auto tau_fac_dnr = -0.78 / pow(dvars.lambdar[il], 2) * lambdar_dnr +
                             0.32 * pow(schmidt_num, 1. / 3.) * sqrt(1. / visc_over_rho) * v_evap_dlambda * lambdar_dnr;
    const auto tau_fac_dqr = -0.78 / pow(dvars.lambdar[il], 2) * lambdar_dqr +
                             0.32 * pow(schmidt_num, 1. / 3.) * sqrt(1. / visc_over_rho) * v_evap_dlambda * lambdar_dqr;

    const auto inv_tau = 2. * constants.pi * state.nr[il] * dvars.rho_dry[il] * dv * tau_fac;
    const auto tau = 1.0 / inv_tau;

    const auto tau_dT = -pow(tau, 2) * (2. * constants.pi * state.nr[il]) *
                        (rho_dry_dT * dv * tau_fac +
                         dvars.rho_dry[il] * dv_dT * tau_fac +
                         dvars.rho_dry[il] * dv * tau_fac_dT);
    const auto tau_dq = -pow(tau, 2) * (2. * constants.pi * state.nr[il] * dv *
                                        (rho_dry_dq * tau_fac + dvars.rho_dry[il] * tau_fac_dq));
    const auto tau_nr = -pow(tau, 2) * (2. * constants.pi * dvars.rho_dry[il] * dv *
                                        (tau_fac + state.nr[il] * tau_fac_dnr));
    const auto tau_qr = -pow(tau, 2) * (2. * constants.pi * dvars.rho_dry[il] * state.nr[il] * dv * tau_fac_dqr);

    const auto dqdt = (abl * tau * q_sat_dry_dT - (q_sat_dry - state.q[il]) * (abl * tau_dT + abl_dT * tau)) / pow(abl * tau, 2);
    const auto dqdq = (-tau - (q_sat_dry - state.q[il]) * tau_dq) * pow(inv_tau, 2) / abl;
    const auto dqdnr = dvars.lambdar[il] == 0.0 ? 0 : (-(q_sat_dry - state.q[il]) / abl * pow(inv_tau, 2) * tau_nr);
    const auto dqdqr = -(q_sat_dry - state.q[il]) / abl * pow(inv_tau, 2) * tau_qr;

    const auto i_t = il;
    const auto i_q = i_t + grid.nlev;
    const auto i_nr = i_q + grid.nlev;
    const auto i_qr = i_nr + grid.nlev;

    elem(i_t, i_t) -= dqdt * constants.latvap / constants.cp;
    elem(i_t, i_q) -= dqdq * constants.latvap / constants.cp;
    elem(i_t, i_nr) -= dqdnr * constants.latvap / constants.cp;
    elem(i_t, i_qr) -= dqdqr * constants.latvap / constants.cp;

    elem(i_q, i_t) += dqdt;
    elem(i_q, i_q) += dqdq;
    elem(i_q, i_nr) += dqdnr;
    elem(i_q, i_qr) += dqdqr;

    elem(i_nr, i_t) -= dqdt * state.nr[il] / state.qr[il];
    elem(i_nr, i_q) -= dqdq * state.nr[il] / state.qr[il];
    elem(i_nr, i_nr) -= (q_sat_dry - state.q[il]) * inv_tau / (abl * state.qr[il]) + (state.nr[il] / state.qr[il]) * dqdnr;
    elem(i_nr, i_qr) -= -(q_sat_dry - state.q[il]) * inv_tau * state.nr[il] / (abl * pow(state.qr[il], 2)) + (state.nr[il] / state.qr[il]) * dqdqr;

    elem(i_qr, i_t) -= dqdt;
    elem(i_qr, i_q) -= dqdq;
    elem(i_qr, i_nr) -= dqdnr;
    elem(i_qr, i_qr) -= dqdqr;
  }
}

double Evaporation::calc_v_evap(const RainshaftConstants &constants, double lambdar) const
{
  if (v_table.has_value())
  {
    double d_micron = 1.e6 / lambdar;
    return v_table->lookup_value(d_micron);
  }
  else
  {
    if (use_numerical_integration)
    {
      return calc_v_evap_numerical(constants, lambdar);
    }
    else
    {
      return calc_v_evap_gamma(constants, lambdar);
    }
  }
}

double Evaporation::calc_v_evap_dlambda(const RainshaftConstants& constants, double lambdar) const {
  if (dv_table.has_value()) {
    double d_micron = 1.e6 / lambdar;
    return dv_table->lookup_value(d_micron);
  } else {
    // if (use_numerical_integration) {
    //   return calc_v_evap_numerical(constants, lambdar);
    // } else {
    //   return calc_v_evap_gamma(constants, lambdar);
    // }

    // currently no numerical integration method
    return calc_v_evap_dlambda_gamma(constants, lambdar);
  }
}

double Evaporation::calc_v_evap_dlambda_gamma(const RainshaftConstants& constants, double lambdar) const {
  // Skip function when no rain present.
  if (lambdar == 0.) {
    return 0.;
  }
  double v_evapi(0.), v_evapi1(0.);
  // Factor converting D^3 to drop mass in grams.
  double d3_to_gram = 1000. * constants.pi * constants.rhow / 6.;
  // Gram conversion factor to various powers.
  double d2g_2third = pow(d3_to_gram, 2./3.);
  double d2g_1third = pow(d3_to_gram, 1./3.);
  double d2g_1sixth = pow(d3_to_gram, 1./6.);
  // Integral for D <= 134.43 micron.
  v_evapi = sqrt(4579.5 * d2g_2third)
    * tgamma_lower(3.5, lambdar * 1.3443e-4) * pow(lambdar, -2.5);
  // Integral for 134.43 micron < D <= 1511.64 micron.
  v_evapi += sqrt(49.62 * d2g_1third)
    * (tgamma(3., lambdar * 1.3443e-4) - tgamma(3., lambdar * 1.51164e-3))
    / (lambdar * lambdar);
  // Integral for 1511.64 micron < D <= 3477.84 micron.
  v_evapi += sqrt(17.32 * d2g_1sixth)
    * (tgamma(2.75, lambdar * 1.51164e-3) - tgamma(2.75, lambdar * 3.47784e-3))
    * pow(lambdar, -1.75);
  // Integral for 3477.84 micron < D.
  v_evapi += sqrt(9.17) * tgamma(2.5, lambdar * 3.47784e-3) * pow(lambdar, -1.5);


  // Integral for D <= 134.43 micron.
  v_evapi1 = sqrt(4579.5 * d2g_2third)
    * tgamma_lower(4.5, lambdar * 1.3443e-4) * pow(lambdar, -3.5);
  // Integral for 134.43 micron < D <= 1511.64 micron.
  v_evapi1 += sqrt(49.62 * d2g_1third)
    * (tgamma(4., lambdar * 1.3443e-4) - tgamma(4., lambdar * 1.51164e-3))
    * pow(lambdar, -3);
  // Integral for 1511.64 micron < D <= 3477.84 micron.
  v_evapi1 += sqrt(17.32 * d2g_1sixth)
    * (tgamma(3.75, lambdar * 1.51164e-3) - tgamma(3.75, lambdar * 3.47784e-3))
    * pow(lambdar, -2.75);
  // Integral for 3477.84 micron < D.
  v_evapi1 += sqrt(9.17) * tgamma(3.5, lambdar * 3.47784e-3) * pow(lambdar, -2.5);
  return v_evapi / lambdar - v_evapi1;
}
