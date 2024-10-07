#ifndef EVAPORATION_HPP
#define EVAPORATION_HPP
#include <optional>
#include <boost/math/special_functions/gamma.hpp>
using boost::math::tgamma, boost::math::tgamma_lower;

#include "lookup_linear.hpp"
#include "rainshaft_process.hpp"
#include "rainshaft_types.hpp"
#include "saturation.hpp"

class Evaporation : public RainshaftProcess
{
private:
  template <bool WithGrad = false>
  RealOptGrad<WithGrad, 1> calc_diffusivity(const double t, const double p) const {
    const double fac = 8.794e-5 / p;
    const double diffusivity = fac * pow(t, 1.81);

    if constexpr (WithGrad) {
      // Derivative with respect to t
      return {diffusivity, {1.81 * diffusivity / t}};
    } else {
      return diffusivity;
    }
  }

  template <bool WithGrad = false>
  RealOptGrad<WithGrad, 1> calc_dynamic_viscosity(const double t) const {
    const double viscosity_over_t = 1.496e-6 * sqrt(t) / (t + 120.);
    const double viscosity = viscosity_over_t * t;

    if constexpr (WithGrad) {
      return {
        viscosity,
        // Derivative with respect to t
        {0.5 * viscosity_over_t * (t + 360.) / (t + 120.)}
      };
    } else {
      return viscosity;
    }
  }

  template <bool WithGrad = false>
  RealOptGrad<WithGrad, 2> calc_visc_over_rho(const double t, const RealOptGrad<WithGrad, 2> rho_dry) const {
    const RealOptGrad<WithGrad, 1> visc = calc_dynamic_viscosity<WithGrad>(t);
    const double visc_over_rho = get_val(visc) / get_val(rho_dry);

    if constexpr (WithGrad) {
      const auto [visc_dT] = get_grad(visc);
      const auto [rho_dry_dT, rho_dry_dq] = get_grad(rho_dry);
      return {visc_over_rho, {
        // Derivative with respect to t
        (visc_dT - visc_over_rho * rho_dry_dT) / get_val(rho_dry),
        // Derivative with respect to q
        -visc_over_rho * rho_dry_dq / get_val(rho_dry)
      }};
    } else {
      return visc_over_rho;
    }
  }

  template <bool WithGrad = false>
  RealOptGrad<WithGrad, 2> calc_schmidt_num(const RealOptGrad<WithGrad, 1> diffusivity, const RealOptGrad<WithGrad, 2> visc_over_rho) const {
    const double schmidt_num = get_val(visc_over_rho) / get_val(diffusivity);

    if constexpr (WithGrad) {
      const auto [dv_dT] = get_grad(diffusivity);
      const auto [visc_over_rho_dT, visc_over_rho_dq] = get_grad(visc_over_rho);
      return {schmidt_num, {
        // Derivative with respect to t
        (visc_over_rho_dT - get_val(visc_over_rho) * dv_dT / get_val(diffusivity)) / get_val(diffusivity),
        // Derivative with respect to q
        visc_over_rho_dq / get_val(diffusivity)
      }};
    } else {
      return schmidt_num;
    }
  }

  template <bool WithGrad = false>
  RealOptGrad<WithGrad, 1> calc_abl(const RainshaftConstants &constants, const double t, const RealOptGrad<WithGrad, 1> q_sat_dry) const {
    const double fac = 1.0 / (constants.cp * constants.rvapor);
    const double l_over_t_2 = pow(constants.latvap / t, 2);
    const double abl = 1. + fac * l_over_t_2 * get_val(q_sat_dry);

    if constexpr (WithGrad) {
      const auto [q_sat_dry_dT] = get_grad(q_sat_dry);
      // Derivative with respect to t
      return {abl, {fac * l_over_t_2 * (q_sat_dry_dT - 2.0 * get_val(q_sat_dry) / t)}};
    } else {
      return abl;
    }
  }

  template <bool WithGrad = false>
  RealOptGrad<WithGrad, 4> calc_tau_inv(const RainshaftConstants &constants,
                            const double nr,
                            const RealOptGrad<WithGrad, 1> diffusivity,
                            const RealOptGrad<WithGrad, 2> rho_dry,
                            const RealOptGrad<WithGrad, 2> visc_over_rho,
                            const RealOptGrad<WithGrad, 2> schmidt_num,
                            const RealOptGrad<WithGrad, 1> v_evap,
                            const RealOptGrad<WithGrad, 2> lambdar) const {
    const double two_pi = 2 * constants.pi;
    const double two_pi_nr = two_pi * nr;
    const double rho_diffusivity = get_val(rho_dry) * get_val(diffusivity);
    const double t1 = two_pi_nr * rho_diffusivity;

    const double scale_inv_lambdar = 0.78 / get_val(lambdar);
    const double sqrt_visc_over_rho = sqrt(get_val(visc_over_rho));
    const double t2_term = 0.32 * cbrt(get_val(schmidt_num)) / sqrt_visc_over_rho;
    const double scaled_evap = t2_term * get_val(v_evap);
    const double t2 = scale_inv_lambdar + scaled_evap;
    
    const double tau_inv = t1 * t2;

    if constexpr (WithGrad) {
      const auto [diffusivity_dT] = get_grad(diffusivity);
      const auto [rho_dry_dT, rho_dry_dq] = get_grad(rho_dry);
      const double t1_dT = two_pi_nr * (get_val(rho_dry) * diffusivity_dT + rho_dry_dT * get_val(diffusivity));
      const double t1_dq = two_pi_nr * get_val(diffusivity) * rho_dry_dq;
      const double t1_dnr = two_pi * rho_diffusivity;

      const auto [schmidt_num_dT, schmidt_num_dq] = get_grad(schmidt_num);
      const auto [visc_over_rho_dT, visc_over_rho_dq] = get_grad(visc_over_rho);
      const double t2_dT = scaled_evap * (schmidt_num_dT / (3.0 * get_val(schmidt_num)) - visc_over_rho_dT / (2.0 * get_val(visc_over_rho)));
      const double t2_dq = scaled_evap * (schmidt_num_dq / (3.0 * get_val(schmidt_num)) - visc_over_rho_dq / (2.0 * get_val(visc_over_rho)));

      const auto [lambdar_dnr, lambdar_dqr] = get_grad(lambdar);
      const auto [v_evap_lambdar] = get_grad(v_evap);
      const double scale_lambdar_neg_two = -scale_inv_lambdar / get_val(lambdar);
      const double t2_dnr = (scale_lambdar_neg_two + t2_term * v_evap_lambdar) * lambdar_dnr;
      const double t2_dqr = (scale_lambdar_neg_two + t2_term * v_evap_lambdar) * lambdar_dqr;

      return {tau_inv, {
        t1 * t2_dT + t2 * t1_dT, // Derivative with respect to t
        t1 * t2_dq + t2 * t1_dq, // Derivative with respect to q
        t1 * t2_dnr + t2 * t1_dnr, // Derivative with respect to nr
        t1 * t2_dqr // Derivative with respect to qr
      }};
    } else {
      return tau_inv;
    }
  }

  template <bool WithGrad = false>
  RealOptGrad<WithGrad, 4> calc_q_evap(const double q, const RealOptGrad<WithGrad, 1> q_sat_dry, const RealOptGrad<WithGrad, 1> abl, const RealOptGrad<WithGrad, 4> tau_inv) const {
    const double subsat = get_val(q_sat_dry) - q;
    const double q_evap_num = subsat * get_val(tau_inv);
    const double q_evap = q_evap_num / get_val(abl);

    if constexpr (WithGrad) {
      const auto [tau_inv_dT, tau_inv_dq, tau_inv_dnr, tau_inv_dqr] = get_grad(tau_inv);
      const auto [abl_dT] = get_grad(abl);
      const auto [q_sat_dry_dT] = get_grad(q_sat_dry);
      return {q_evap, {
        (subsat * (tau_inv_dT - get_val(tau_inv) * abl_dT / get_val(abl)) + get_val(tau_inv) * q_sat_dry_dT) / get_val(abl),
        (subsat * tau_inv_dq - get_val(tau_inv)) / get_val(abl),
        subsat * tau_inv_dnr / get_val(abl),
        subsat * tau_inv_dqr / get_val(abl)
      }};
    } else {
      return q_evap;
    }
  }

  template <bool WithGrad = false>
  RealOptGrad<WithGrad, 4> calc_n_evap(const double nr, const double qr, const RealOptGrad<WithGrad, 4> q_evap) const
  {
    const double nr_over_qr = nr / qr;
    const double n_evap = nr_over_qr * get_val(q_evap);

    if constexpr (WithGrad) {
      const auto [q_evap_dT, q_evap_dq, q_evap_dnr, q_evap_dqr] = get_grad(q_evap);
      return {n_evap, {
        nr_over_qr * q_evap_dT,
        nr_over_qr * q_evap_dq,
        nr_over_qr * q_evap_dnr + get_val(q_evap) / qr,
        nr_over_qr * (q_evap_dqr - get_val(q_evap) / qr)
      }};
    } else {
      return n_evap;
    }
  }

  template <bool WithGrad = false>
  RealOptGrad<WithGrad, 4> calc_T_evap(const RainshaftConstants &constants, const RealOptGrad<WithGrad, 4> q_evap) const
  {
    const double fac = -constants.latvap / constants.cp;
    const double t_tend = fac * get_val(q_evap);

    if constexpr (WithGrad) {
      const auto [q_evap_dT, q_evap_dq, q_evap_dnr, q_evap_dqr] = get_grad(q_evap);
      return {t_tend, {
        fac * q_evap_dT,
        fac * q_evap_dq,
        fac * q_evap_dnr,
        fac * q_evap_dqr,
      }};
    } else {
      return t_tend;
    }
  }

public:
  // SPS: Should use a different pointer type to guarantee this
  // does not outlast the sat_form object.
  Evaporation(const RainshaftConstants &constants,
              const SaturationFormulae &sat_form_in, bool use_v_table,
              bool use_numerical_integration);

  // Calculate tendency from current state.
  void calc_tend(const RainshaftConstants &constants,
                              const RainshaftGrid &grid,
                              const StateConst& state,
                              const RainshaftDerivedVars &dvars,
                              Tendency& tend) const;

  void calc_tend_jac_prod(const RainshaftConstants &constants,
                          const RainshaftGrid &grid,
                          const StateConst& state,
                          const RainshaftDerivedVars &dvars,
                          const double *const vec,
                          double *const prod) const;

  void calc_tend_jac(const RainshaftConstants &constants,
                             const RainshaftGrid &grid,
                             const StateConst& state,
                             const RainshaftDerivedVars &dvars,
                             Matrix jac) const;

  // Calculate characteristic velocity used for velocity calculation.
  template <bool WithGrad = false>
  RealOptGrad<WithGrad, 1> calc_v_evap(const RainshaftConstants &constants, double lambdar) const
  {
    if (v_table.has_value())
    {
      const double d_micron = 1.e6 / lambdar;
      const RealGrad<1> v_evap = v_table->lookup_value(d_micron);
      if constexpr (WithGrad) {
        return {get_val(v_evap), {-get_grad(v_evap)[0] * d_micron / lambdar}};
      } else {
        return get_val(v_evap);
      }
    }
    else
    {
      return calc_v_evap<WithGrad>(constants, lambdar, use_numerical_integration);
    }
  }

  template <bool WithGrad = false>
std::array<RealOptGrad<WithGrad, 4>, 3> calc_evap(const RainshaftConstants& constants, const double t, const double q,
  const double nr, const double qr, const double p_mid, const RealOptGrad<WithGrad, 2> rho_dry, const RealOptGrad<WithGrad, 2> lambdar) const
  {
    // Skip the rest of this if no rain.
    if (qr < constants.qsmall)
    {
      return {};
    }
    const RealOptGrad<WithGrad, 1> q_sat_dry = sat_form.q_sat_dry<WithGrad>(t, p_mid);
    // Skip the rest of this if not saturated.
    if (get_val(q_sat_dry) < q)
    {
      return {};
    }

    const RealOptGrad<WithGrad, 1> diffusivity = calc_diffusivity<WithGrad>(t, p_mid);
    const RealOptGrad<WithGrad, 2> visc_over_rho = calc_visc_over_rho<WithGrad>(t, rho_dry);
    const RealOptGrad<WithGrad, 2> schmidt_num = calc_schmidt_num<WithGrad>(diffusivity, visc_over_rho);
    const RealOptGrad<WithGrad, 1> abl = calc_abl<WithGrad>(constants, t, q_sat_dry);
    const RealOptGrad<WithGrad, 1> v_evap = calc_v_evap<WithGrad>(constants, get_val(lambdar));
    const RealOptGrad<WithGrad, 4> tau_inv = calc_tau_inv<WithGrad>(constants, nr, diffusivity, rho_dry, visc_over_rho, schmidt_num, v_evap, lambdar);
    const RealOptGrad<WithGrad, 4> q_evap = calc_q_evap<WithGrad>(q, q_sat_dry, abl, tau_inv);
    const RealOptGrad<WithGrad, 4> n_evap = calc_n_evap<WithGrad>(nr, qr, q_evap);
    const RealOptGrad<WithGrad, 4> t_evap = calc_T_evap<WithGrad>(constants, q_evap);

    return {q_evap, n_evap, t_evap};
  }

  const bool use_numerical_integration;

private:
  // Water vapor saturation formulae.
  const SaturationFormulae &sat_form;

  // Particle velocity lookup table.
  // This is currently hard-coded to P3 settings, i.e. it contains 20 entries
  // for values every 10 microns between 5 and 195 micron, and 280 entries for
  // values every 30 microns between 195 and 8595 micron.
  const std::optional<LookupLinear> v_table;

  static std::optional<LookupLinear> create_lookup(const RainshaftConstants &constants,
                                      const bool use_v_table,
                                      const bool use_numerical_integration);

  template <bool WithGrad = false>
  static RealOptGrad<WithGrad, 1> calc_v_evap(const RainshaftConstants &constants, const double lambdar, const bool use_numerical_integration)
  {
    if (use_numerical_integration) {
      return calc_v_evap_numerical<WithGrad>(constants, lambdar);
    } else {
      return calc_v_evap_gamma<WithGrad>(constants, lambdar);
    }
  }

  // Calculate characteristic velocity used for velocity calculation.
  // This version ignores the lookup table, if present, and always just
  // calculates using incomplete gamma functions.
  template <bool WithGrad = false>
  static RealOptGrad<WithGrad, 1> calc_v_evap_gamma(const RainshaftConstants &constants, const double lambdar)
  {
    // Skip function when no rain present.
    if (lambdar == 0.)
    {
      // This returns 0 (and 0 derivative)
      return {};
    }
    // Factor converting D^3 to drop mass in grams.
    const double d3_to_gram = 1000. * constants.pi * constants.rhow / 6.;
    
    // Integral for D <= 134.43 micron.
    const double int1_fac = sqrt(4579.5) * cbrt(d3_to_gram);
    const double lambdar_neg_2_5 = pow(lambdar, -2.5);
    const double lambdar_size1 = lambdar * 1.3443e-4;
    const double term1 = int1_fac * tgamma_lower(3.5, lambdar_size1) * lambdar_neg_2_5;

    // Integral for 134.43 micron < D <= 1511.64 micron.
    const double int2_fac = sqrt(49.62) * pow(d3_to_gram, 1.0 / 6.0);
    const double lambdar_neg_2 = pow(lambdar, -2);
    const double lambdar_size2 = lambdar * 1.51164e-3;
    const double term2 = int2_fac * (tgamma(3., lambdar_size1) - tgamma(3., lambdar_size2)) * lambdar_neg_2;

    // Integral for 1511.64 micron < D <= 3477.84 micron.
    const double int3_fac = sqrt(17.32) * pow(d3_to_gram, 1.0 / 12.0);
    const double lambdar_neg_1_75 = pow(lambdar, -1.75);
    const double lambdar_size3 = lambdar * 3.47784e-3;
    const double term3 = int3_fac * (tgamma(2.75, lambdar_size2) - tgamma(2.75, lambdar_size3)) * lambdar_neg_1_75;
    
    // Integral for 3477.84 micron < D.
    const double int4_fac = sqrt(9.17);
    const double lambdar_neg_1_5 = pow(lambdar, -1.5);
    const double term4 = int4_fac * tgamma(2.5, lambdar_size3) * lambdar_neg_1_5;

    const double v_evap = term1 + term2 + term3 + term4;

    if constexpr (WithGrad) {
      const double term1_dl = (term1 - int1_fac * tgamma_lower(4.5, lambdar_size1) * lambdar_neg_2_5) / lambdar;
      const double term2_dl = (term2 - int2_fac * (tgamma(4., lambdar_size1) - tgamma(4., lambdar_size2)) * lambdar_neg_2) / lambdar;
      const double term3_dl = (term3 - int3_fac * (tgamma(3.75, lambdar_size2) - tgamma(3.75, lambdar_size3)) * lambdar_neg_1_75) / lambdar;
      const double term4_dl = (term4 - int4_fac * tgamma(3.5, lambdar_size3) * lambdar_neg_1_5) / lambdar;
      return {v_evap, {term1_dl + term2_dl + term3_dl + term4_dl}};
    } else {
      return v_evap;
    }
  }

  // Calculate characteristic velocity (v_evap) using Riemann sum over midpoints
  template <bool WithGrad = false>
  static RealOptGrad<WithGrad, 1> calc_v_evap_numerical(const RainshaftConstants &constants, const double lambdar)
  {
    constexpr double dd = 2.; // micron
    constexpr std::size_t low_k = 1;
    constexpr std::size_t high_k = 10000;
    const double amg_fac = 1000. * constants.pi * constants.rhow / 6.;
    RealOptGrad<WithGrad, 1> accum{};
    for (std::size_t kk = low_k; kk != high_k + 1; ++kk) {
      const double dia_micron = (kk - 0.5) * dd;
      const double dia = dia_micron * 1.e-6;
      const double vt = [=] () {
        const double amg = amg_fac * pow(dia, 3);
        if (dia_micron <= 134.43) {
          return 4.5795e3 * pow(amg, 2./3.);
        } else if (dia_micron <= 1511.64) {
          return 4.962e1 * pow(amg, 1./3.);
        } else if (dia_micron <= 3477.84) {
          return 1.732e1 * pow(amg, 1./6.);
        } else {
          return 9.17;
        }
      }();
      const double integrand = sqrt(vt * dia) * dia * exp(-lambdar * dia);
      get_val(accum) += integrand;

      if constexpr (WithGrad) {
        get_grad(accum)[0] -= dia * integrand;
      }
    }

    // Multiply by quadrature cell width and unit conversion
    const double integral_scale = dd * 1.e-6;
    if constexpr (WithGrad) {
      // Apply product rule to lambdar * int_0^inf f(D, lambdar) dD
      get_grad(accum)[0] = integral_scale * (lambdar * get_grad(accum)[0] + get_val(accum));
    }
    // Multiply by lambdar factored outside integral
    get_val(accum) *= lambdar * integral_scale;
    return accum;
  }
};

#endif // EVAPORATION_HPP
