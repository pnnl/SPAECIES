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

  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const auto [q_evap, n_evap, t_evap] = calc_evap(constants, t[il], q[il], nr[il], qr[il], grid.p_mid[il], dvars.rho_dry[il], dvars.lambdar[il]);
    t_tend[il] = t_evap;
    q_tend[il] = q_evap;
    nr_tend[il] = -n_evap;
    qr_tend[il] = -q_evap;
  }
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
    const RealGrad<2> rho_dry = dvars.get_rho_dry<true>(constants, t[il], q[il], il);
    const RealGrad<2> lambdar = dvars.get_lambdar<true>(constants, nr[il], qr[il], il);
    const auto [q_evap, n_evap, t_evap] = calc_evap<true>(constants, t[il], q[il], nr[il], qr[il], grid.p_mid[il], rho_dry, lambdar);

    const std::size_t i_t = il;
    const std::size_t i_q = i_t + grid.nlev;
    const std::size_t i_nr = i_q + grid.nlev;
    const std::size_t i_qr = i_nr + grid.nlev;

    const auto [t_evap_dT, t_evap_dq, t_evap_dnr, t_evap_dqr] = get_grad(t_evap);
    jac(i_t, i_t) += t_evap_dT;
    jac(i_t, i_q) += t_evap_dq;
    jac(i_t, i_nr) += t_evap_dnr;
    jac(i_t, i_qr) += t_evap_dqr;

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
