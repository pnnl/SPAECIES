#include "evaporation.hpp"
#include "saturation.hpp"
#include <cstddef>
#include <cmath>

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
    : use_numerical_integration(use_numerical_integration), sat_form(sat_form_in),
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

    const auto [n_evap_dT, n_evap_dq, n_evap_dnr, n_evap_dqr] = get_grad(n_evap);
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
