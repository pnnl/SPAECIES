#include "evaporation.hpp"
#include "saturation.hpp"
#include <cstddef>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
using boost::math::tgamma, boost::math::tgamma_lower;
using std::sqrt, std::pow;

// Calculate v_evap using incomplete gamma functions.
double Evaporation::calc_v_evap_gamma(const RainshaftConstants &constants, double lambdar)
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
double Evaporation::calc_v_evap_numerical(const RainshaftConstants &constants, double lambdar)
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

double Evaporation::calc_v_evap(const RainshaftConstants &constants, const double lambdar, const bool use_numerical_integration)
{
  if (use_numerical_integration) {
    return calc_v_evap_numerical(constants, lambdar);
  } else {
    return calc_v_evap_gamma(constants, lambdar);
  }
}

std::optional<LookupLinear> Evaporation::create_lookup(const RainshaftConstants &constants,
                                      const bool use_v_table,
                                      const bool use_numerical_integration)
{
  if (use_v_table) {
    return std::make_optional<LookupLinear>(
      std::vector{5., 195., 8595.},
      std::vector{1., 1.},
      [=] (const double x) {
        return calc_v_evap(constants, 1.e6 / x, use_numerical_integration);
      });
  } else {
    return std::nullopt;
  }
}

Evaporation::Evaporation(const RainshaftConstants &constants,
                         const SaturationFormulae &sat_form_in,
                         const bool use_v_table,
                         const bool use_numerical_integration)
    : sat_form(sat_form_in), use_numerical_integration(use_numerical_integration),
    v_table(create_lookup(constants, use_v_table, use_numerical_integration))
{
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
    const auto q_sat_dry = sat_form.q_sat_dry(state.t[il], grid.p_mid[il]);
    // Skip the rest of this if not saturated.
    if (q_sat_dry < state.q[il])
    {
      continue;
    }

    const auto diffusivity = calc_diffusivity(state.t[il], grid.p_mid[il]);
    const auto rho_dry = dvars.get_rho_dry(constants, state, il);
    const auto visc_over_rho = calc_visc_over_rho(state.t[il], rho_dry);
    const auto schmidt_num = calc_schmidt_num(diffusivity, visc_over_rho);
    const auto abl = calc_abl(constants, state.t[il], q_sat_dry);
    const auto lambdar = dvars.get_lambdar(constants, state, il);
    const auto v_evap = calc_v_evap(constants, dvars.lambdar[il]);
    const auto tau_inv = calc_tau_inv(constants, state.nr[il], diffusivity, rho_dry, visc_over_rho, schmidt_num, v_evap, lambdar);
    const auto q_evap = calc_q_evap(state.q[il], q_sat_dry, abl, tau_inv);
    
    t_tend[il] = calc_T_tend(constants, q_evap);
    q_tend[il] = q_evap;
    nr_tend[il] = -calc_n_evap(state.nr[il], state.qr[il], q_evap);
    qr_tend[il] = -q_evap;
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
                                Matrix jac) const
{
  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    // Skip the rest of this if no rain.
    if (state.qr[il] < constants.qsmall)
    {
      continue;
    }
    const auto q_sat_dry = sat_form.q_sat_dry<true>(state.t[il], grid.p_mid[il]);

    // Skip the rest of this if not saturated.
    if (get_val(q_sat_dry) < state.q[il])
    {
      continue;
    }

    const auto diffusivity = calc_diffusivity<true>(state.t[il], grid.p_mid[il]);
    const auto rho_dry = dvars.get_rho_dry<true>(constants, state, il);
    const auto visc_over_rho = calc_visc_over_rho<true>(state.t[il], rho_dry);
    const auto schmidt_num = calc_schmidt_num<true>(diffusivity, visc_over_rho);
    const auto abl = calc_abl<true>(constants, state.t[il], q_sat_dry);
    const auto lambdar = dvars.get_lambdar<true>(constants, state, il);
    const ValGrad<1> v_evap = {calc_v_evap(constants, dvars.lambdar[il]), {0.0}}; // TODO: fix deriv
    const auto tau_inv = calc_tau_inv<true>(constants, state.nr[il], diffusivity, rho_dry, visc_over_rho, schmidt_num, v_evap, lambdar);
    const auto q_evap = calc_q_evap<true>(state.q[il], q_sat_dry, abl, tau_inv);
    const auto n_evap = calc_n_evap<true>(state.nr[il], state.qr[il], q_evap);
    const auto t_tend = calc_T_tend<true>(constants, q_evap);

    const auto i_t = il;
    const auto i_q = i_t + grid.nlev;
    const auto i_nr = i_q + grid.nlev;
    const auto i_qr = i_nr + grid.nlev;

    const auto [t_tend_dT, t_tend_dq, t_tend_dnr, t_tend_dqr] = get_grad(t_tend);
    jac(i_t, i_t) += t_tend_dT;
    jac(i_t, i_q) += t_tend_dq;
    jac(i_t, i_nr) += t_tend_dnr;
    jac(i_t, i_qr) += t_tend_dqr;

    const auto [q_evap_dT, q_evap_dq, q_evap_dnr, q_evap_dqr] = get_grad(q_evap);
    jac(i_q, i_t) += q_evap_dT;
    jac(i_q, i_q) += q_evap_dq;
    jac(i_q, i_nr) += q_evap_dnr;
    jac(i_q, i_qr) += q_evap_dqr;

    const auto [n_evap_dT, n_evap_dq, n_evap_dnr, n_evap_dqr] = get_grad(q_evap);
    jac(i_nr, i_t) -= n_evap_dT;
    jac(i_nr, i_q) -= n_evap_dq;
    jac(i_nr, i_nr) -= n_evap_dnr;
    jac(i_nr, i_qr) -= n_evap_dqr;

    jac(i_qr, i_t) -= q_evap_dT;
    jac(i_qr, i_q) -= q_evap_dq;
    jac(i_qr, i_nr) -= q_evap_dnr;
    jac(i_qr, i_qr) -= q_evap_dqr;
  }
}

double Evaporation::calc_v_evap(const RainshaftConstants &constants, double lambdar) const
{
  if (v_table.has_value())
  {
    double d_micron = 1.e6 / lambdar;
    return get_val(v_table->lookup_value(d_micron));
  }
  else
  {
    return calc_v_evap(constants, lambdar, use_numerical_integration);
  }
}
