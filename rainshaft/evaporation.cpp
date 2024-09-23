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
      [=] (const auto x) {
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

void Evaporation::calc_tend(const RainshaftConstants& constants,
                            const RainshaftGrid& grid,
                            const StateConst& state,
                            const RainshaftDerivedVars& dvars,
                            Tendency& tend) const {
  VarConst t = state.get_variable("T");
  VarConst q = state.get_variable("q");
  VarConst nr = state.get_variable("nr");
  VarConst qr = state.get_variable("qr");
  VarMut t_tend = tend.get_variable("T_tend");
  VarMut q_tend = tend.get_variable("q_tend");
  VarMut nr_tend = tend.get_variable("nr_tend");
  VarMut qr_tend = tend.get_variable("qr_tend");

  for (std::size_t il = 0; il != grid.nlev; ++il) {
    // Skip the rest of this if no rain.
    if (qr[il] < constants.qsmall) {
      continue;
    }
    const auto q_sat_dry = sat_form.q_sat_dry(t[il], grid.p_mid[il]);
    // Skip the rest of this if not saturated.
    if (q_sat_dry < q[il]) {
      continue;
    }

    const auto diffusivity = calc_diffusivity(t[il], grid.p_mid[il]);
    const auto rho_dry = dvars.get_rho_dry(constants, state, il);
    const auto visc_over_rho = calc_visc_over_rho(t[il], rho_dry);
    const auto schmidt_num = calc_schmidt_num(diffusivity, visc_over_rho);
    const auto abl = calc_abl(constants, t[il], q_sat_dry);
    const auto lambdar = dvars.get_lambdar(constants, state, il);
    const auto v_evap = calc_v_evap(constants, dvars.lambdar[il]);
    const auto tau_inv = calc_tau_inv(constants, nr[il], diffusivity, rho_dry, visc_over_rho, schmidt_num, v_evap, lambdar);
    const auto q_evap = calc_q_evap(q[il], q_sat_dry, abl, tau_inv);
    
    t_tend[il] = calc_T_tend(constants, q_evap);
    q_tend[il] = q_evap;
    nr_tend[il] = -calc_n_evap(nr[il], qr[il], q_evap);
    qr_tend[il] = -q_evap;
  }
}

void Evaporation::calc_tend_jac_prod(const RainshaftConstants &constants,
                                     const RainshaftGrid &grid,
                                     const StateConst& state,
                                     const RainshaftDerivedVars &dvars,
                                     const double *const vec,
                                     double *const prod) const
{
  throw std::logic_error("Jacobian product not implemented");
}

void Evaporation::calc_tend_jac(const RainshaftConstants &constants,
                                const RainshaftGrid &grid,
                                const StateConst& state,
                                const RainshaftDerivedVars &dvars,
                                Matrix jac) const
{
  VarConst t = state.get_variable("T");
  VarConst q = state.get_variable("q");
  VarConst nr = state.get_variable("nr");
  VarConst qr = state.get_variable("qr");

  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    // Skip the rest of this if no rain.
    if (qr[il] < constants.qsmall)
    {
      continue;
    }
    const auto q_sat_dry = sat_form.q_sat_dry<true>(t[il], grid.p_mid[il]);
    // Skip the rest of this if not saturated.
    if (get_val(q_sat_dry) < q[il])
    {
      continue;
    }

    const auto diffusivity = calc_diffusivity<true>(t[il], grid.p_mid[il]);
    const auto rho_dry = dvars.get_rho_dry<true>(constants, state, il);
    const auto visc_over_rho = calc_visc_over_rho<true>(t[il], rho_dry);
    const auto schmidt_num = calc_schmidt_num<true>(diffusivity, visc_over_rho);
    const auto abl = calc_abl<true>(constants, t[il], q_sat_dry);
    const auto lambdar = dvars.get_lambdar<true>(constants, state, il);
    const auto v_evap = calc_v_evap<true>(constants, dvars.lambdar[il]);    
    const auto tau_inv = calc_tau_inv<true>(constants, nr[il], diffusivity, rho_dry, visc_over_rho, schmidt_num, v_evap, lambdar);
    const auto q_evap = calc_q_evap<true>(q[il], q_sat_dry, abl, tau_inv);
    const auto n_evap = calc_n_evap<true>(nr[il], qr[il], q_evap);
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
