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
      std::vector{5., 8595., 360000.},
      std::vector{0.001, 1.},
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

double Sedimentation::calc_max_step(const RainshaftConstants &constants,
                                    const RainshaftGrid &grid,
                                    const RainshaftDerivedVars& dvars,
                                    double cfl) const
{
  double lambdar_top = std::cbrt(constants.pi * constants.rhow
                                  * constants.nr_top / constants.qr_top);
  const auto [v0, v3] = rain_fall_speeds(constants, dvars.rho_dry[0],
                                      lambdar_top);
  double spatial_frequency = v3 / dvars.dz[0];
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    const auto [v0, v3] = rain_fall_speeds(constants, dvars.rho_dry[il],
                                    dvars.lambdar[il]);
    spatial_frequency = std::max(spatial_frequency, v3 / dvars.dz[il]);
  }

  return cfl / spatial_frequency;
}

void Sedimentation::calc_tend(const RainshaftConstants &constants,
                              const RainshaftGrid &grid,
                              const StateConst &state,
                              const RainshaftDerivedVars &dvars,
                              Tendency& tend) const
{
  const double lambdar_top = cbrt(constants.pi * constants.rhow * constants.nr_top / constants.qr_top);
  auto [v0_prev, v3_prev] = rain_fall_speeds(constants, constants.rho_top, lambdar_top);
  double nr_prev = constants.nr_top;
  double qr_prev = constants.qr_top;
  double rho_prev = constants.rho_top;
  double nr_flux_prev = calc_nr_flux(nr_prev, rho_prev, v0_prev);
  double qr_flux_prev = calc_qr_flux(qr_prev, rho_prev, v3_prev);

  VarConst nr = *state.get_variable("nr");
  VarConst qr = *state.get_variable("qr");
  VarMut nr_tend = *tend.get_variable("nr_tend");
  VarMut qr_tend = *tend.get_variable("qr_tend");

  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const auto [v0, v3] = rain_fall_speeds(constants, dvars.rho_dry[il], dvars.lambdar[il]);
    const double nr_flux = calc_nr_flux(nr[il], dvars.rho_dry[il], v0);
    const double qr_flux = calc_qr_flux(qr[il], dvars.rho_dry[il], v3);
    
    nr_tend[il] += calc_nr_tend(dvars.dz[il], dvars.rho_dry[il], nr_flux, nr_flux_prev);
    qr_tend[il] += calc_qr_tend(dvars.dz[il], dvars.rho_dry[il], qr_flux, qr_flux_prev);

    rho_prev = dvars.rho_dry[il];
    nr_flux_prev = nr_flux;
    qr_flux_prev = qr_flux;
  }
}

void Sedimentation::calc_tend_jac(const RainshaftConstants &constants,
                                  const RainshaftGrid &grid,
                                  const StateConst& state,
                                  const RainshaftDerivedVars &dvars,
                                  Matrix jac) const
{
  const double lambdar_top = cbrt(constants.pi * constants.rhow * constants.nr_top / constants.qr_top);
  auto [v0_prev, v3_prev] = rain_fall_speeds<true>(constants, {constants.rho_top, {0., 0.}}, {lambdar_top, {0., 0.}});
  double nr_prev = constants.nr_top;
  double qr_prev = constants.qr_top;
  RealGrad<2> rho_prev{constants.rho_top, {0., 0.}};
  RealGrad<4> nr_flux_prev = calc_nr_flux<true>(nr_prev, rho_prev, v0_prev);
  RealGrad<4> qr_flux_prev = calc_qr_flux<true>(qr_prev, rho_prev, v3_prev);

  VarConst t = *state.get_variable("T");
  VarConst q = *state.get_variable("q");
  VarConst nr = *state.get_variable("nr");
  VarConst qr = *state.get_variable("qr");

  const std::size_t offset_t = *state.get_idx("T");
  const std::size_t offset_q = *state.get_idx("q");
  const std::size_t offset_nr = *state.get_idx("nr");
  const std::size_t offset_qr = *state.get_idx("qr");

  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const RealGrad<2> rho = dvars.get_rho_dry<true>(constants, t[il], q[il], il);
    const RealGrad<2> lambdar = dvars.get_lambdar<true>(constants, nr[il], qr[il], il);
    const RealGrad<2> dz = dvars.get_dz<true>(constants, t[il], q[il], il);
    const auto [v0, v3] = rain_fall_speeds<true>(constants, rho, lambdar);
    const RealGrad<4> nr_flux = calc_nr_flux<true>(nr[il], rho, v0);
    const RealGrad<4> qr_flux = calc_qr_flux<true>(qr[il], rho, v3);

    const RealGrad<8> nr_tend = calc_nr_tend<true>(dz, rho, nr_flux, nr_flux_prev);
    const RealGrad<8> qr_tend = calc_qr_tend<true>(dz, rho, qr_flux, qr_flux_prev);
    const auto [nr_tend_dT_prev, nr_tend_dT, nr_tend_dq_prev, nr_tend_dq, nr_tend_dnr_prev, nr_tend_dnr, nr_tend_dqr_prev, nr_tend_dqr] = get_grad(nr_tend);
    const auto [qr_tend_dT_prev, qr_tend_dT, qr_tend_dq_prev, qr_tend_dq, qr_tend_dnr_prev, qr_tend_dnr, qr_tend_dqr_prev, qr_tend_dqr] = get_grad(qr_tend);

    const std::size_t i_t = offset_t + il;
    const std::size_t i_q = offset_q + il;
    const std::size_t i_nr = offset_nr + il;
    const std::size_t i_qr = offset_qr + il;

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

std::set<std::string> Sedimentation::get_required_vars() const
{
  return {"T", "q", "nr", "qr"};
}

std::set<std::string> Sedimentation::get_optional_vars() const
{
  return {};
}
