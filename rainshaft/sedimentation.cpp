#include "sedimentation.hpp"
#include <cstddef>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
using boost::math::tgamma, boost::math::tgamma_lower;
using std::pow, std::sqrt, std::cbrt, std::exp;
std::optional<LookupLinear> Sedimentation::create_lookup(const RainshaftConstants &constants,
                                      const bool use_v_table,
                                      const bool create_v0,
                                      const bool use_numerical_integration)
{
  if (use_v_table) {
    return std::make_optional<LookupLinear>(
      std::vector{5., 195., 8595.},
      std::vector{1., 1.},
      [=] (const double x) {
          const auto [v0, v3] = rain_fall_speeds(constants, 1.e6 / x, use_numerical_integration);
          return create_v0 ? v0 : v3;
      });
  } else {
    return std::nullopt;
  }
}

Sedimentation::Sedimentation(const RainshaftConstants &constants, bool use_v_table,
                             bool use_numerical_integration)
    : use_numerical_integration(use_numerical_integration), 
      v0_table(create_lookup(constants, use_v_table, true, use_numerical_integration)),
      v3_table(create_lookup(constants, use_v_table, false, use_numerical_integration))
{}

void Sedimentation::calc_tend(const RainshaftConstants &constants,
                              const RainshaftGrid &grid,
                              const StateConst &state,
                              const RainshaftDerivedVars &dvars,
                              Tendency& tend) const
{
  const auto lambdar_top = cbrt(constants.pi * constants.rhow * constants.nr_top / constants.qr_top);
  auto [v0_prev, v3_prev] = rain_fall_speeds(constants, constants.rho_top, lambdar_top);
  auto nr_prev = constants.nr_top;
  auto qr_prev = constants.qr_top;
  auto rho_prev = constants.rho_top;
  auto nr_flux_prev = calc_nr_flux(nr_prev, rho_prev, v0_prev);
  auto qr_flux_prev = calc_qr_flux(qr_prev, rho_prev, v3_prev);

  VarConst nr = state.get_variable("nr");
  VarConst qr = state.get_variable("qr");
  VarMut nr_tend = tend.get_variable("nr_tend");
  VarMut qr_tend = tend.get_variable("qr_tend");

  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const auto [v0, v3] = rain_fall_speeds(constants, dvars.rho_dry[il], dvars.lambdar[il]);
    const auto nr_flux = calc_nr_flux(nr[il], dvars.rho_dry[il], v0);
    const auto qr_flux = calc_qr_flux(qr[il], dvars.rho_dry[il], v3);
    
    nr_tend[il] = calc_nr_tend(dvars.dz[il], dvars.rho_dry[il], rho_prev, nr_flux, nr_flux_prev);
    qr_tend[il] = calc_qr_tend(dvars.dz[il], dvars.rho_dry[il], rho_prev, qr_flux, qr_flux_prev);

    rho_prev = dvars.rho_dry[il];
    nr_flux_prev = nr_flux;
    qr_flux_prev = qr_flux;
  }
}

void Sedimentation::calc_tend_jac_prod(const RainshaftConstants &constants,
                                       const RainshaftGrid &grid,
                                       const StateConst& state,
                                       const RainshaftDerivedVars &dvars,
                                       const double *const vec,
                                       double *const prod) const
{
  throw std::logic_error("Jacobian product not implemented");
}

void Sedimentation::calc_tend_jac(const RainshaftConstants &constants,
                                  const RainshaftGrid &grid,
                                  const StateConst& state,
                                  const RainshaftDerivedVars &dvars,
                                  Matrix jac) const
{
  const auto lambdar_top = cbrt(constants.pi * constants.rhow * constants.nr_top / constants.qr_top);
  auto [v0_prev, v3_prev] = rain_fall_speeds<true>(constants, {constants.rho_top, {0., 0.}}, {lambdar_top, {0., 0.}});
  auto nr_prev = constants.nr_top;
  auto qr_prev = constants.qr_top;
  RealOptGrad<true, 2> rho_prev{constants.rho_top, {0., 0.}};
  auto nr_flux_prev = calc_nr_flux<true>(nr_prev, rho_prev, v0_prev);
  auto qr_flux_prev = calc_qr_flux<true>(qr_prev, rho_prev, v3_prev);

  VarConst t = state.get_variable("T");
  VarConst q = state.get_variable("q");
  VarConst nr = state.get_variable("nr");
  VarConst qr = state.get_variable("qr");

  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const auto rho = dvars.get_rho_dry<true>(constants, t[il], q[il], il);
    const auto lambdar = dvars.get_lambdar<true>(constants, nr[il], qr[il], il);
    const auto dz = dvars.get_dz<true>(constants, t[il], q[il], il);
    const auto [v0, v3] = rain_fall_speeds<true>(constants, rho, lambdar);
    const auto nr_flux = calc_nr_flux<true>(nr[il], rho, v0);
    const auto qr_flux = calc_qr_flux<true>(qr[il], rho, v3);

    const auto nr_tend = calc_nr_tend<true>(dz, rho, rho_prev, nr_flux, nr_flux_prev);
    const auto qr_tend = calc_qr_tend<true>(dz, rho, rho_prev, qr_flux, qr_flux_prev);
    const auto [nr_tend_dT_prev, nr_tend_dT, nr_tend_dq_prev, nr_tend_dq, nr_tend_dnr_prev, nr_tend_dnr, nr_tend_dqr_prev, nr_tend_dqr] = get_grad(nr_tend);
    const auto [qr_tend_dT_prev, qr_tend_dT, qr_tend_dq_prev, qr_tend_dq, qr_tend_dnr_prev, qr_tend_dnr, qr_tend_dqr_prev, qr_tend_dqr] = get_grad(qr_tend);

    const auto i_t = il;
    const auto i_q = i_t + grid.nlev;
    const auto i_nr = i_q + grid.nlev;
    const auto i_qr = i_nr + grid.nlev;

    jac(i_nr, i_t) += nr_tend_dT;
    jac(i_nr, i_q) += nr_tend_dq;
    jac(i_nr, i_nr) += nr_tend_dnr;
    jac(i_nr, i_qr) += nr_tend_dqr;

    jac(i_qr, i_t) += qr_tend_dT;
    jac(i_qr, i_q) += qr_tend_dq;
    jac(i_qr, i_nr) += qr_tend_dnr;
    jac(i_qr, i_qr) += qr_tend_dqr;

    if (il > 0) {
      jac(i_nr, i_t - 1) += nr_tend_dT_prev;
      jac(i_nr, i_q - 1) += nr_tend_dq_prev;
      jac(i_nr, i_nr - 1) += nr_tend_dnr_prev;
      jac(i_nr, i_qr - 1) += nr_tend_dqr_prev;

      jac(i_qr, i_t - 1) += qr_tend_dT_prev;
      jac(i_qr, i_q - 1) += qr_tend_dq_prev;
      jac(i_qr, i_nr - 1) += qr_tend_dnr_prev;
      jac(i_qr, i_qr - 1) += qr_tend_dqr_prev;
    }

    rho_prev = rho;
    nr_flux_prev = nr_flux;
    qr_flux_prev = qr_flux;
  }
}
