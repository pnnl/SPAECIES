#ifndef EVAPORATION_HPP
#define EVAPORATION_HPP
#include <vector>
#include <optional>
#include <boost/math/special_functions/gamma.hpp>
using boost::math::tgamma, boost::math::tgamma_lower;

#include "lookup_linear.hpp"
#include "rainshaft_process.hpp"
#include "saturation.hpp"

class Evaporation : public RainshaftProcess
{
private:
  template <bool WithGrad = false>
  Val<WithGrad, 1> calc_diffusivity(const double t, const double p) const {
    const auto fac = 8.794e-5 / p;
    const auto diffusivity = fac * pow(t, 1.81);

    if constexpr (WithGrad) {
      // Derivative with respect to t
      return {diffusivity, {1.81 * fac * pow(t, 0.81)}};
    } else {
      return diffusivity;
    }
  }

  template <bool WithGrad = false>
  Val<WithGrad, 1> calc_dynamic_viscosity(const double t) const {
    const auto fac = 1.496e-6 / (t + 120.);
    const auto viscosity = fac * pow(t, 1.5);

    if constexpr (WithGrad) {
      return {
        viscosity,
        // Derivative with respect to t
        {0.5 * fac * sqrt(t) * (t + 360.) / (t + 120.)}
      };
    } else {
      return viscosity;
    }
  }

  template <bool WithGrad = false>
  Val<WithGrad, 2> calc_visc_over_rho(const double t, const Val<WithGrad, 2> rho_dry) const {
    const auto visc = calc_dynamic_viscosity<WithGrad>(t);
    const auto visc_over_rho = get_val(visc) / get_val(rho_dry);

    if constexpr (WithGrad) {
      const auto [visc_dT] = get_grad(visc);
      const auto [rho_dry_dT, rho_dry_dq] = get_grad(rho_dry);
      return {visc_over_rho, {
        // Derivative with respect to t
        (visc_dT - get_val(visc) * rho_dry_dT / get_val(rho_dry)) / get_val(rho_dry),
        // Derivative with respect to q
         -get_val(visc) * pow(get_val(rho_dry), -2) * rho_dry_dq
      }};
    } else {
      return visc_over_rho;
    }
  }

  template <bool WithGrad = false>
  Val<WithGrad, 2> calc_schmidt_num(const Val<WithGrad, 1> diffusivity, const Val<WithGrad, 2> visc_over_rho) const {
    const auto schmidt_num = get_val(visc_over_rho) / get_val(diffusivity);

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
  Val<WithGrad, 1> calc_abl(const RainshaftConstants &constants, const double t, const Val<WithGrad, 1> q_sat_dry) const {
    const auto fac = 1.0 / (constants.cp * constants.rvapor);
    const auto l_over_t_2 = pow(constants.latvap / t, 2);
    const auto abl = 1. + fac * l_over_t_2 * get_val(q_sat_dry);

    if constexpr (WithGrad) {
      const auto [q_sat_dry_dT] = get_grad(q_sat_dry);
      // Derivative with respect to t
      return {abl, {fac * l_over_t_2 * (q_sat_dry_dT - 2.0 * get_val(q_sat_dry) / t)}};
    } else {
      return abl;
    }
  }

  template <bool WithGrad = false>
  Val<WithGrad, 4> calc_tau_inv(const RainshaftConstants &constants,
                            const double nr,
                            const Val<WithGrad, 1> diffusivity,
                            const Val<WithGrad, 2> rho_dry,
                            const Val<WithGrad, 2> visc_over_rho,
                            const Val<WithGrad, 2> schmidt_num,
                            const Val<WithGrad, 1> v_evap,
                            const Val<WithGrad, 2> lambdar) const {
    const auto two_pi = 2 * constants.pi;
    const auto two_pi_nr = two_pi * nr;
    const auto rho_diffusivity = get_val(rho_dry) * get_val(diffusivity);
    const auto t1 = two_pi_nr * rho_diffusivity;

    const auto scale_inv_lambdar = 0.78 / get_val(lambdar);
    const auto sqrt_visc_over_rho = sqrt(get_val(visc_over_rho));
    const auto t2_term = 0.32 * cbrt(get_val(schmidt_num)) / sqrt_visc_over_rho;
    const auto scaled_evap = t2_term * get_val(v_evap);
    const auto t2 = scale_inv_lambdar + scaled_evap;
    
    const auto tau_inv = t1 * t2;

    if constexpr (WithGrad) {
      const auto [diffusivity_dT] = get_grad(diffusivity);
      const auto [rho_dry_dT, rho_dry_dq] = get_grad(rho_dry);
      const auto t1_dT = two_pi_nr * (get_val(rho_dry) * diffusivity_dT + rho_dry_dT * get_val(diffusivity));
      const auto t1_dq = two_pi_nr * get_val(diffusivity) * rho_dry_dq;
      const auto t1_dnr = two_pi * rho_diffusivity;

      const auto [schmidt_num_dT, schmidt_num_dq] = get_grad(schmidt_num);
      const auto [visc_over_rho_dT, visc_over_rho_dq] = get_grad(visc_over_rho);
      const auto t2_dT = scaled_evap * (schmidt_num_dT / (3.0 * get_val(schmidt_num)) - visc_over_rho_dT / (2.0 * get_val(visc_over_rho)));
      const auto t2_dq = scaled_evap * (schmidt_num_dq / (3.0 * get_val(schmidt_num)) - visc_over_rho_dq / (2.0 * get_val(visc_over_rho)));

      const auto [lambdar_dnr, lambdar_dqr] = get_grad(lambdar);
      const auto [v_evap_lambdar] = get_grad(v_evap);
      const auto scale_lambdar_neg_two = -scale_inv_lambdar / get_val(lambdar);
      const auto t2_dnr = (scale_lambdar_neg_two + t2_term * v_evap_lambdar) * lambdar_dnr;
      const auto t2_dqr = (scale_lambdar_neg_two + t2_term * v_evap_lambdar) * lambdar_dqr;

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
  Val<WithGrad, 4> calc_q_evap(const double q, const Val<WithGrad, 1> q_sat_dry, const Val<WithGrad, 1> abl, const Val<WithGrad, 4> tau_inv) const {
    const auto q_diff = get_val(q_sat_dry) - q;
    const auto q_evap_num = q_diff * get_val(tau_inv);
    const auto q_evap = q_evap_num / get_val(abl);

    if constexpr (WithGrad) {
      const auto [tau_inv_dT, tau_inv_dq, tau_inv_dnr, tau_inv_dqr] = get_grad(tau_inv);
      const auto [abl_dT] = get_grad(abl);
      const auto [q_sat_dry_dT] = get_grad(q_sat_dry);
      return {q_evap, {
        (q_diff * (tau_inv_dT - get_val(tau_inv) * abl_dT / get_val(abl)) + get_val(tau_inv) * q_sat_dry_dT) / get_val(abl),
        (q_diff * tau_inv_dq - get_val(tau_inv)) / get_val(abl),
        q_diff * tau_inv_dnr / get_val(abl),
        q_diff * tau_inv_dqr / get_val(abl)
      }};
    } else {
      return q_evap;
    }
  }

  template <bool WithGrad = false>
  Val<WithGrad, 4> calc_n_evap(const double nr, const double qr, const Val<WithGrad, 4> q_evap) const
  {
    const auto qr_over_nr = nr / qr;
    const auto n_evap = qr_over_nr * get_val(q_evap);

    if constexpr (WithGrad) {
      const auto [q_evap_dT, q_evap_dq, q_evap_dnr, q_evap_dqr] = get_grad(q_evap);
      return {n_evap, {
        qr_over_nr * q_evap_dT,
        qr_over_nr * q_evap_dq,
        qr_over_nr * q_evap_dnr + get_val(q_evap) / qr,
        qr_over_nr * (q_evap_dqr - get_val(q_evap) / qr)
      }};
    } else {
      return n_evap;
    }
  }

  template <bool WithGrad = false>
  Val<WithGrad, 4> calc_T_tend(const RainshaftConstants &constants, const Val<WithGrad, 4> q_evap) const
  {
    const auto fac = -constants.latvap / constants.cp;
    const auto t_tend = fac * get_val(q_evap);

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
  RainshaftTendency calc_tend(const RainshaftConstants &constants,
                              const RainshaftGrid &grid,
                              const RainshaftState &state,
                              const RainshaftDerivedVars &dvars) const;

  void calc_tend_jac_prod(const RainshaftConstants &constants,
                          const RainshaftGrid &grid,
                          const RainshaftState &state,
                          const RainshaftDerivedVars &dvars,
                          const double *const vec,
                          double *const prod) const;

  void calc_tend_jac(const RainshaftConstants &constants,
                             const RainshaftGrid &grid,
                             const RainshaftState &state,
                             const RainshaftDerivedVars &dvars,
                             Matrix jac) const;

  // Calculate characteristic velocity used for velocity calculation.
  template <bool WithGrad = false>
  Val<WithGrad, 1> calc_v_evap(const RainshaftConstants &constants, double lambdar) const
  {
    if (v_table.has_value())
    {
      const auto d_micron = 1.e6 / lambdar;
      const auto v_evap = v_table->lookup_value(d_micron);
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
  static Val<WithGrad, 1> calc_v_evap(const RainshaftConstants &constants, const double lambdar, const bool use_numerical_integration)
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
  static Val<WithGrad, 1> calc_v_evap_gamma(const RainshaftConstants &constants, const double lambdar)
  {
    // Skip function when no rain present.
    if (lambdar == 0.)
    {
      // This returns 0 (and 0 derivative)
      return {};
    }
    // Factor converting D^3 to drop mass in grams.
    const auto d3_to_gram = 1000. * constants.pi * constants.rhow / 6.;
    
    // Integral for D <= 134.43 micron.
    const auto int1_fac = sqrt(4579.5) * cbrt(d3_to_gram);
    const auto lambdar_neg_2_5 = pow(lambdar, -2.5);
    const auto lambdar_size1 = lambdar * 1.3443e-4;
    const auto term1 = int1_fac * tgamma_lower(3.5, lambdar_size1) * lambdar_neg_2_5;

    // Integral for 134.43 micron < D <= 1511.64 micron.
    const auto int2_fac = sqrt(49.62) * pow(d3_to_gram, 1.0 / 6.0);
    const auto lambdar_neg_2 = pow(lambdar, -2);
    const auto lambdar_size2 = lambdar * 1.51164e-3;
    const auto term2 = int2_fac * (tgamma(3., lambdar_size1) - tgamma(3., lambdar_size2)) * lambdar_neg_2;

    // Integral for 1511.64 micron < D <= 3477.84 micron.
    const auto int3_fac = sqrt(17.32) * pow(d3_to_gram, 1.0 / 12.0);
    const auto lambdar_neg_1_75 = pow(lambdar, -1.75);
    const auto lambdar_size3 = lambdar * 3.47784e-3;
    const auto term3 = int3_fac * (tgamma(2.75, lambdar_size2) - tgamma(2.75, lambdar_size3)) * lambdar_neg_1_75;
    
    // Integral for 3477.84 micron < D.
    const auto int4_fac = sqrt(9.17);
    const auto lambdar_neg_1_5 = pow(lambdar, -1.5);
    const auto term4 = int4_fac * tgamma(2.5, lambdar_size3) * lambdar_neg_1_5;

    const auto v_evap = term1 + term2 + term3 + term4;

    if constexpr (WithGrad) {
      const auto term1_dl = (term1 - int1_fac * tgamma_lower(4.5, lambdar_size1) * lambdar_neg_2_5) / lambdar;
      const auto term2_dl = (term2 - int2_fac * (tgamma(4., lambdar_size1) - tgamma(4., lambdar_size2)) * lambdar_neg_2) / lambdar;
      const auto term3_dl = (term3 - int3_fac * (tgamma(3.75, lambdar_size2) - tgamma(3.75, lambdar_size3)) * lambdar_neg_1_75) / lambdar;
      const auto term4_dl = (term4 - int4_fac * tgamma(3.5, lambdar_size3) * lambdar_neg_1_5) / lambdar;
      return {v_evap, {term1_dl + term2_dl + term3_dl + term4_dl}};
    } else {
      return v_evap;
    }
  }

  // Calculate characteristic velocity (v_evap) using Riemann sum over midpoints
  template <bool WithGrad = false>
  static Val<WithGrad, 1> calc_v_evap_numerical(const RainshaftConstants &constants, const double lambdar)
  {
    constexpr auto dd = 2.; // micron
    constexpr std::size_t low_k = 1;
    constexpr std::size_t high_k = 10000;
    const auto amg_fac = 1000. * constants.pi * constants.rhow / 6.;
    Val<WithGrad, 1> accum{};
    for (std::size_t kk = low_k; kk != high_k + 1; ++kk) {
      const auto dia_micron = (kk - 0.5) * dd;
      const auto dia = dia_micron * 1.e-6;
      const auto vt = [=] () {
        const auto amg = amg_fac * pow(dia, 3);
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
      get_val(accum) += sqrt(vt * dia) * dia * exp(-lambdar * dia);

      if constexpr (WithGrad) {
        get_grad(accum)[0] -= dia * get_val(accum);
      }
    }

    // Multiply by quadrature cell width, unit conversion, and lambdar factored out of integral
    const auto scale = dd * 1.e-6 * lambdar;
    get_val(accum) *= scale;
    if constexpr (WithGrad) {
      // Apply product rule to lambdar * int_0^inf f(D, lambdar) dD
      get_grad(accum)[0] = scale * get_grad(accum)[0] + get_val(accum);
    }
    return accum;
  }
};

#endif // EVAPORATION_HPP
