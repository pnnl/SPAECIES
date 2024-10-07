#ifndef SEDIMENTATION_HPP
#define SEDIMENTATION_HPP
#include <vector>
#include <optional>
#include <tuple>
#include <iostream>
#include <boost/math/special_functions/gamma.hpp>
using boost::math::tgamma, boost::math::tgamma_lower;

#include "lookup_linear.hpp"
#include "rainshaft_process.hpp"

class Sedimentation : public RainshaftProcess 
{
private:
  static std::optional<LookupLinear> create_lookup(const RainshaftConstants &constants,
                                      const bool use_v_table,
                                      const bool create_v0,
                                      const bool use_numerical_integration);

  template<bool WithGrad = false>
  RealOptGrad<WithGrad, 4> calc_nr_flux(const double nr,
                                const RealOptGrad<WithGrad, 2> rho,
                                const RealOptGrad<WithGrad, 4> v0) const
  {
    const double nr_flux = get_val(v0) * nr * get_val(rho);

    if constexpr (WithGrad) {
      const double nr_times_rho = nr * get_val(rho);
      const double v0_times_nr = get_val(v0) * nr;

      const double nr_flux_dT = get_grad(v0)[0] * nr_times_rho + v0_times_nr * get_grad(rho)[0];
      const double nr_flux_dq = get_grad(v0)[1] * nr_times_rho + v0_times_nr * get_grad(rho)[1];
      const double nr_flux_dnr = get_grad(v0)[2] * nr_times_rho + get_val(v0) * get_val(rho);
      const double nr_flux_dqr = get_grad(v0)[3] * nr_times_rho;

      return {nr_flux, {nr_flux_dT, nr_flux_dq, nr_flux_dnr, nr_flux_dqr}};
    } else {
      return nr_flux;
    }
  }

  template<bool WithGrad = false>
  RealOptGrad<WithGrad, 4> calc_qr_flux(const double qr,
                                const RealOptGrad<WithGrad, 2> rho,
                                const RealOptGrad<WithGrad, 4> v3) const
  {
    const double qr_flux = get_val(v3) * qr * get_val(rho);

    if constexpr (WithGrad) {
      const double qr_times_rho = qr * get_val(rho);
      const double v3_times_qr = get_val(v3) * qr;

      const double qr_flux_dT = get_grad(v3)[0] * qr_times_rho + v3_times_qr * get_grad(rho)[0];
      const double qr_flux_dq = get_grad(v3)[1] * qr_times_rho + v3_times_qr * get_grad(rho)[1];
      const double qr_flux_dnr = get_grad(v3)[2] * qr_times_rho;
      const double qr_flux_dqr = get_grad(v3)[3] * qr_times_rho + get_val(v3) * get_val(rho);

      return {qr_flux, {qr_flux_dT, qr_flux_dq, qr_flux_dnr, qr_flux_dqr}};
    } else {
      return qr_flux;
    }
  }

  template<bool WithGrad = false>
  RealOptGrad<WithGrad, 8> calc_nr_tend(const RealOptGrad<WithGrad, 2> dz,
                                const RealOptGrad<WithGrad, 2> rho,
                                const RealOptGrad<WithGrad, 2> rho_prev,
                                const RealOptGrad<WithGrad, 4> nr_flux_in,
                                const RealOptGrad<WithGrad, 4> nr_flux_prev_in) const
  {
    const double nr_flux_prev = get_val(nr_flux_prev_in);
    const double nr_flux = get_val(nr_flux_in);
    const double dz_times_rho = get_val(dz) * get_val(rho);

    const double nr_tend = (nr_flux_prev - nr_flux) / dz_times_rho;

    if constexpr (WithGrad) {
      const auto [nr_flux_prev_dT, nr_flux_prev_dq, nr_flux_prev_dnr, nr_flux_prev_dqr] = get_grad(nr_flux_prev_in);
      const auto [nr_flux_dT, nr_flux_dq, nr_flux_dnr, nr_flux_dqr] = get_grad(nr_flux_in);

      const double tend_dTprev = nr_flux_prev_dT / dz_times_rho;
      const double tend_dT = (-nr_flux_dT - (nr_flux_prev - nr_flux) * (get_val(dz) * get_grad(rho)[0] + get_grad(dz)[0] * get_val(rho)) / dz_times_rho) / dz_times_rho;

      const double tend_dqprev = nr_flux_prev_dq / dz_times_rho;
      const double tend_dq = (-nr_flux_dq - (nr_flux_prev - nr_flux) * (get_val(dz) * get_grad(rho)[1] + get_grad(dz)[1] * get_val(rho)) / dz_times_rho) / dz_times_rho;

      const double tend_dnrprev = nr_flux_prev_dnr / dz_times_rho;
      const double tend_dnr = -nr_flux_dnr / dz_times_rho;

      const double tend_dqrprev = nr_flux_prev_dqr / dz_times_rho;
      const double tend_dqr = -nr_flux_dqr / dz_times_rho;

      return {nr_tend, {tend_dTprev, tend_dT, tend_dqprev, tend_dq, tend_dnrprev, tend_dnr, tend_dqrprev, tend_dqr}};
    } else {
      return nr_tend;
    }
  }

  template<bool WithGrad = false>
  RealOptGrad<WithGrad, 8> calc_qr_tend(const RealOptGrad<WithGrad, 2> dz,
                                const RealOptGrad<WithGrad, 2> rho,
                                const RealOptGrad<WithGrad, 2> rho_prev,
                                const RealOptGrad<WithGrad, 4> qr_flux_in,
                                const RealOptGrad<WithGrad, 4> qr_flux_prev_in) const
  {
    const double qr_flux_prev = get_val(qr_flux_prev_in);
    const double qr_flux = get_val(qr_flux_in);
    const double dz_times_rho = get_val(dz) * get_val(rho);

    const double qr_tend = (qr_flux_prev - qr_flux) / dz_times_rho;

    if constexpr (WithGrad) {
      const auto [qr_flux_prev_dT, qr_flux_prev_dq, qr_flux_prev_dnr, qr_flux_prev_dqr] = get_grad(qr_flux_prev_in);
      const auto [qr_flux_dT, qr_flux_dq, qr_flux_dnr, qr_flux_dqr] = get_grad(qr_flux_in);

      const double tend_dTprev = qr_flux_prev_dT / dz_times_rho;
      const double tend_dT = (-qr_flux_dT - (qr_flux_prev - qr_flux) * (get_val(dz) * get_grad(rho)[0] + get_grad(dz)[0] * get_val(rho)) / dz_times_rho) / dz_times_rho;

      const double tend_dqprev = qr_flux_prev_dq / dz_times_rho;
      const double tend_dq = (-qr_flux_dq - (qr_flux_prev - qr_flux) * (get_val(dz) * get_grad(rho)[1] + get_grad(dz)[1] * get_val(rho)) / dz_times_rho) / dz_times_rho;

      const double tend_dnrprev = qr_flux_prev_dnr / dz_times_rho;
      const double tend_dnr = -qr_flux_dnr / dz_times_rho;

      const double tend_dqrprev = qr_flux_prev_dqr / dz_times_rho;
      const double tend_dqr = -qr_flux_dqr / dz_times_rho;

      return {qr_tend, {tend_dTprev, tend_dT, tend_dqprev, tend_dq, tend_dnrprev, tend_dnr, tend_dqrprev, tend_dqr}};
    } else {
      return qr_tend;
    }
  }

public:
  template<bool WithGrad = false, int I = 4>
  using Speeds = std::tuple<RealOptGrad<WithGrad, I>, RealOptGrad<WithGrad, I>>;

  // Constructor allowing use of lookup table.
  Sedimentation(const RainshaftConstants& constants, bool use_v_table,
                bool use_numerical_integration);

  // Calculate tendency from current state.
  void calc_tend(const RainshaftConstants& constants,
                 const RainshaftGrid& grid,
                 const StateConst& state,
                 const RainshaftDerivedVars& dvars,
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

  // For a given value of lambdar, what are the rain number and mass fall speeds
  // and their derivatives wrt (T, q, nr, qr)?
  template <bool WithGrad = false>
  Speeds<WithGrad, 4> rain_fall_speeds(const RainshaftConstants &constants, 
                                           RealOptGrad<WithGrad, 2> rho, 
                                           RealOptGrad<WithGrad, 2> lambdar) const
  {
    // Calculate correction for air density, with reference state at 1000 hPa and 273.15 K.
    const RealOptGrad<WithGrad, 2> rf = rho_fac<WithGrad>(constants, rho);

    if (v0_table.has_value() && v3_table.has_value())
    {
      // skip computation if no rain is present
      if (get_val(lambdar) == 0.)
      {
        // This returns 0 (and 0 derivative)
        return {};
      }

      // Extract values from lookup table. Note: derivatives wrt lambdar are computed
      // from finite differences in lookup_value, NOT 
      const double d_micron = 1.e6 / get_val(lambdar);
      const RealGrad<1> v0 = v0_table->lookup_value(d_micron); 
      const RealGrad<1> v3 = v3_table->lookup_value(d_micron);

      if constexpr (WithGrad) {
        return {{get_val(rf) * get_val(v0), {get_grad(rf)[0] * get_val(v0),
                                             get_grad(rf)[1] * get_val(v0), 
                                             get_val(rf) * -get_grad(v0)[0] * d_micron / get_val(lambdar) * get_grad(lambdar)[0],
                                             get_val(rf) * -get_grad(v0)[0] * d_micron / get_val(lambdar) * get_grad(lambdar)[1]}},
                {get_val(rf) * get_val(v3), {get_grad(rf)[0] * get_val(v3), 
                                             get_grad(rf)[1] * get_val(v3), 
                                             get_val(rf) * -get_grad(v3)[0] * d_micron / get_val(lambdar) * get_grad(lambdar)[0],
                                             get_val(rf) * -get_grad(v3)[0] * d_micron / get_val(lambdar) * get_grad(lambdar)[1]}}};
      } else {
        return {get_val(rf) * get_val(v0), get_val(rf) * get_val(v3)};
      }
    }
    else
    {
      const auto [v0, v3] = rain_fall_speeds<WithGrad>(constants, get_val(lambdar), use_numerical_integration);

      if constexpr (WithGrad) {
        return {{get_val(rf) * get_val(v0), {get_grad(rf)[0] * get_val(v0),
                                             get_grad(rf)[1] * get_val(v0), 
                                             get_val(rf) * get_grad(v0)[0] * get_grad(lambdar)[0],
                                             get_val(rf) * get_grad(v0)[0] * get_grad(lambdar)[1]}},
                {get_val(rf) * get_val(v3), {get_grad(rf)[0] * get_val(v3), 
                                             get_grad(rf)[1] * get_val(v3), 
                                             get_val(rf) * get_grad(v3)[0] * get_grad(lambdar)[0],
                                             get_val(rf) * get_grad(v3)[0] * get_grad(lambdar)[1]}}};
      } else {
        return {get_val(rf) * get_val(v0), get_val(rf) * get_val(v3)};
      }
    }
  }

  const bool use_numerical_integration;

private:

  // Particle velocity lookup tables.
  // This is currently hard-coded to P3 settings, i.e. it contains 20 entries
  // for values every 10 microns between 5 and 195 micron, and 280 entries for
  // values every 30 microns between 195 and 8595 micron.
  std::optional<LookupLinear> v0_table;
  std::optional<LookupLinear> v3_table;

  template <bool WithGrad = false>
  static Speeds<WithGrad, 1> rain_fall_speeds(const RainshaftConstants &constants, 
                                           const double lambdar, 
                                           const bool use_numerical_integration)
  {
    if (use_numerical_integration) {
      return rain_fall_speeds_numerical<WithGrad>(constants, lambdar);
    } else {
      return rain_fall_speeds_gamma<WithGrad>(constants, lambdar);
    }
  }


  // What would the speeds for a given lambdar be at standard temperature and
  // pressure?
  // This version ignores any lookup table, if present, and always returns fall
  // speeds calculated with incomplete gamma functions. Returns the derivative
  // of speeds wrt lambdar only
  template <bool WithGrad = false>
  static Speeds<WithGrad, 1> rain_fall_speeds_gamma(const RainshaftConstants& constants,
                                                 double lambdar)
  {
    // Skip function when no rain present.
    if (lambdar == 0.)
    {
      // This returns 0 (and 0 derivative)
      return {};
    }

    double d3_to_gram = 1000. * constants.pi * constants.rhow / 6.;

    // v0
    const double int1_fac = 4579.5 * pow(d3_to_gram, 2. / 3.);
    const double lambdar_neg_2 = pow(lambdar, -2);
    const double lambdar_size1 = lambdar * 1.3443e-4;
    const double term1v0 = int1_fac * tgamma_lower(3., lambdar_size1) * lambdar_neg_2;

    const double int2_fac = 49.62 * pow(d3_to_gram, 1. / 3.);
    const double lambdar_neg_1 = 1.0 / lambdar;
    const double lambdar_size2 = lambdar * 1.51164e-3;
    const double term2v0 = int2_fac * (tgamma(2., lambdar_size1) - tgamma(2., lambdar_size2)) * lambdar_neg_1;

    const double int3_fac = 17.32 * pow(d3_to_gram, 1. / 6.);
    const double lambdar_neg_0_5 = pow(lambdar, -0.5);
    const double lambdar_size3 = lambdar * 3.47784e-3;
    const double term3v0 = int3_fac * (tgamma(1.5, lambdar_size2) - tgamma(1.5, lambdar_size3)) * lambdar_neg_0_5;

    const double int4_fac = 9.17;
    const double term4v0 = int4_fac * tgamma(1., lambdar_size3);

    const double v0 = term1v0 + term2v0 + term3v0 + term4v0;

    // v3
    const double gamma4 = tgamma(4.);
    const double term1v3 = int1_fac * tgamma_lower(6., lambdar_size1) * lambdar_neg_2 / gamma4;
    const double term2v3 = int2_fac * (tgamma(5., lambdar_size1) - tgamma(5., lambdar_size2)) * lambdar_neg_1 / gamma4;
    const double term3v3 = int3_fac * (tgamma(4.5, lambdar_size2) - tgamma(4.5, lambdar_size3)) * lambdar_neg_0_5 / gamma4;
    const double term4v3 = int4_fac * tgamma(4., lambdar_size3) / gamma4;

    const double v3 = (term1v3 + term2v3 + term3v3 + term4v3);

    if constexpr (WithGrad) {
      const double term1v0_dl = (term1v0 - int1_fac * tgamma_lower(4., lambdar_size1) * lambdar_neg_2) / lambdar;
      const double term2v0_dl = (term2v0 - int2_fac * (tgamma(3., lambdar_size1) - tgamma(3., lambdar_size2)) * lambdar_neg_1) / lambdar;
      const double term3v0_dl = (term3v0 - int3_fac * (tgamma(2.5, lambdar_size2) - tgamma(2.5, lambdar_size3)) * lambdar_neg_0_5) / lambdar;
      const double term4v0_dl = (term4v0 - int4_fac * tgamma(2., lambdar_size3)) / lambdar;

      const double gamma5 = tgamma(5.);
      const double term1v3_dl = (term1v3 - int1_fac * tgamma_lower(7., lambdar_size1) * lambdar_neg_2 / gamma5) * 4.0 / lambdar;
      const double term2v3_dl = (term2v3 - int2_fac * (tgamma(6., lambdar_size1) - tgamma(6., lambdar_size2)) * lambdar_neg_1 / gamma5) * 4.0 / lambdar;
      const double term3v3_dl = (term3v3 - int3_fac * (tgamma(5.5, lambdar_size2) - tgamma(5.5, lambdar_size3)) * lambdar_neg_0_5 / gamma5) * 4.0 / lambdar;
      const double term4v3_dl = (term4v3 - int4_fac * tgamma(5., lambdar_size3) / gamma5) * 4.0 / lambdar;

      return {{v0, {(term1v0_dl + term2v0_dl + term3v0_dl + term4v0_dl)}},
              {v3, {(term1v3_dl + term2v3_dl + term3v3_dl + term4v3_dl)}}};
    } else {
      return {v0, v3};
    }
  }

  // What would the speeds for a given lambdar be at standard temperature and
  // pressure?
  // This version ignores any lookup table, if present, and always returns fall
  // speeds calculated with numerical integration. Returns the derivative
  // of speeds wrt lambdar only
  template <bool WithGrad = false>
  static Speeds<WithGrad, 1> rain_fall_speeds_numerical(const RainshaftConstants& constants,
                                                     double lambdar)
  {
    // Skip function when no rain present.
    if (lambdar == 0.)
    {
      // This returns 0 (and 0 derivative)
      return {};
    }

    constexpr double dd = 2.; // micron
    constexpr std::size_t low_k = 1;
    constexpr std::size_t high_k = 10000;
    const double amg_fac = 1000. * constants.pi * constants.rhow / 6.;
    RealOptGrad<WithGrad, 1> accum0{};
    RealOptGrad<WithGrad, 1> accum3{};
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

      const double integrand0 = vt * exp(-lambdar * dia);
      const double integrand3 = vt * pow(dia, 3) * exp(-lambdar * dia);
      get_val(accum0) += integrand0;
      get_val(accum3) += integrand3;

      if constexpr (WithGrad) {
        get_grad(accum0)[0] -= dia * integrand0;
        get_grad(accum3)[0] -= dia * integrand3;
      } 
    }

    // Multiply by quadrature cell width and unit conversion and scaling factor
    //  int_0^inf exp(-lambdar*D) = 1/lambdar for v0
    //  int_0^inf D^3 * exp(-lambdar*D) = 6/lambdar^4 for v3
    const double integral_scale = dd * 1.e-6;

    if constexpr (WithGrad) {
      // Apply product rule to lambdar * int_0^inf f(D, lambdar) dD for v0 and
      // lambdar^4 / 6 * int_0^inf f(D, lambdar) dD for v3
      get_grad(accum0)[0] = integral_scale * (get_grad(accum0)[0] * lambdar + get_val(accum0));
      get_grad(accum3)[0] = integral_scale * (get_grad(accum3)[0] * pow(lambdar, 4.) / 6.0 + 4.0 * pow(lambdar, 3.) / 6.0 * get_val(accum3));
    }

    get_val(accum0) *= integral_scale * lambdar;
    get_val(accum3) *= integral_scale * pow(lambdar, 4.) / 6.0;

    return {accum0, accum3};
  }


  template <bool WithGrad = false>
  static RealOptGrad<WithGrad, 2> rho_fac(const RainshaftConstants &constants, 
                                  RealOptGrad<WithGrad, 2> rho_dry)
  {
    const double rf = pow(1.e5 / (get_val(rho_dry) * constants.rdry * 273.15), 0.54);
    if constexpr (WithGrad) {
      const double fac = -0.54 * rf / get_val(rho_dry);
      return {rf, {fac * get_grad(rho_dry)[0],
                   fac * get_grad(rho_dry)[1]}};
    } else {
      return rf;
    }
  }

};

#endif // SEDIMENTATION_HPP
