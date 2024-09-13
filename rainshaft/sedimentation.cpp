#include "sedimentation.hpp"
#include <cstddef>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
using std::pow, std::sqrt, std::cbrt, std::exp;
using boost::math::tgamma, boost::math::tgamma_lower;

Sedimentation::Sedimentation(const RainshaftConstants& constants, bool use_v_table,
                             bool use_numerical_integration)
  : use_numerical_integration(use_numerical_integration) {
  if (use_v_table) {
    std::vector<double> range_bounds = {5., 195., 8595.};
    std::vector<double> spacings = {1., 1.};
    std::vector<double> d_microns = LookupTable::calc_x_values(range_bounds,
                                                               spacings);
    std::vector<double> v0_values(d_microns.size(), 0.);
    std::vector<double> v3_values(d_microns.size(), 0.);
    for (std::size_t i = 0; i != d_microns.size(); ++i) {
      double lambdar = 1.e6 / d_microns[i];
      std::tie(v0_values[i], v3_values[i]) = use_numerical_integration ?
        rain_fall_speeds_stp_numerical(constants, lambdar)
        :  rain_fall_speeds_stp_gamma(constants, lambdar);
    }
    v0_table.emplace(range_bounds, spacings, v0_values);
    v3_table.emplace(range_bounds, spacings, v3_values);
  }
}

void Sedimentation::calc_tend(const RainshaftConstants &constants,
                              const RainshaftGrid &grid,
                              const StateConst &state,
                              const RainshaftDerivedVars &dvars,
                              const Tendency& tend) const
{
  const double lambdar_top = cbrt(constants.pi * constants.rhow * constants.nr_top / constants.qr_top);
  const auto [speeds_top0, speeds_top3] = rain_fall_speeds(constants, constants.rho_top, lambdar_top);
  double nr_prev_flux = speeds_top0 * constants.nr_top * constants.rho_top;
  double qr_prev_flux = speeds_top3 * constants.qr_top * constants.rho_top;
  VarConst nr = state.get_variable("nr");
  VarConst qr = state.get_variable("qr");
  VarMut nr_tend = tend.get_variable("nr_tend");
  VarMut qr_tend = tend.get_variable("qr_tend");

  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const auto [speeds0, speeds3] = rain_fall_speeds(constants, dvars.rho_dry[il], dvars.lambdar[il]);
    const double nr_flux = speeds0 * nr[il] * dvars.rho_dry[il];
    const double qr_flux = speeds3 * qr[il] * dvars.rho_dry[il];

    nr_tend[il] = (nr_prev_flux - nr_flux) / (dvars.dz[il] * dvars.rho_dry[il]);
    qr_tend[il] = (qr_prev_flux - qr_flux) / (dvars.dz[il] * dvars.rho_dry[il]);

    nr_prev_flux = nr_flux;
    qr_prev_flux = qr_flux;
  }
}

Sedimentation::Speeds Sedimentation::rain_fall_speeds(const RainshaftConstants& constants,
                                                    double rho_dry, double lambdar) const {
  // Catches case where there is no rain present.
  if (lambdar == 0.) {
    return {0.0, 0.0};
  }
  const auto [v0, v3] = rain_fall_speeds_stp(constants, lambdar);
  // Calculate correction for air density, with reference state at 1000 hPa and 273.15 K.
  const double rf = rho_fac(constants, rho_dry);
  return {v0 * rf, v3 * rf};
}

Sedimentation::Speeds Sedimentation::rain_fall_speeds_stp(const RainshaftConstants& constants,
                                                        double lambdar) const {
  if (v0_table.has_value()) {
    const double d_micron = 1.e6 / lambdar;
    return {v0_table->lookup_value(d_micron), v3_table->lookup_value(d_micron)};
  } else {
    if (use_numerical_integration) {
      return rain_fall_speeds_stp_numerical(constants, lambdar);
    } else {
      return rain_fall_speeds_stp_gamma(constants, lambdar);
    }
  }
}

Sedimentation::Speeds Sedimentation::rain_fall_speeds_stp_gamma(const RainshaftConstants& constants,
                                                              double lambdar) const {
  // Catches case where there is no rain present.
  if (lambdar == 0.) {
    return {0.0, 0.0};
  }
  // Number and mass fall speeds in m/s.
  double v0(0.), v3(0.);
  // Factor converting D^3 to drop mass in grams.
  double d3_to_gram = 1000. * constants.pi * constants.rhow / 6.;
  // Gram conversion factor to various powers.
  double d2g_2third = pow(d3_to_gram, 2./3.);
  double d2g_1third = pow(d3_to_gram, 1./3.);
  double d2g_1sixth = pow(d3_to_gram, 1./6.);
  // Start with number term for D <= 134.43 micron.
  v0 += 4579.5 * d2g_2third * tgamma_lower(3, lambdar * 1.3443e-4) / (lambdar * lambdar);
  // Number term for 134.43 micron < D <= 1511.64 micron.
  v0 += 49.62 * d2g_1third
    * (tgamma(2, lambdar * 1.3443e-4) - tgamma(2, lambdar * 1.51164e-3)) / lambdar;
  // Number term for 1511.64 micron < D <= 3477.84 micron.
  v0 += 17.32 * d2g_1sixth
    * (tgamma(1.5, lambdar * 1.51164e-3) - tgamma(1.5, lambdar * 3.47784e-3)) / sqrt(lambdar);
  // Number term for 3477.84 micron < D.
  v0 += 9.17 * exp(-lambdar * 3.47784e-3);
  // Mass term for D <= 134.43 micron.
  v3 += 4579.5 * d2g_2third * tgamma_lower(6, lambdar * 1.3443e-4) / (lambdar * lambdar);
  // Mass term for 134.43 micron < D <= 1511.64 micron.
  v3 += 49.62 * d2g_1third
    * (tgamma(5, lambdar * 1.3443e-4) - tgamma(5, lambdar * 1.51164e-3)) / lambdar;
  // Mass term for 1511.64 micron < D <= 3477.84 micron.
  v3 += 17.32 * d2g_1sixth
    * (tgamma(4.5, lambdar * 1.51164e-3) - tgamma(4.5, lambdar * 3.47784e-3)) / sqrt(lambdar);
  // Mass term for 3477.84 micron < D.
  v3 += 9.17 * tgamma(4, lambdar * 3.47784e-3);
  // Include the normalization gamma(3) for v3.
  v3 /= 6.;
  return {v0, v3};
}

Sedimentation::Speeds Sedimentation::rain_fall_speeds_stp_numerical(const RainshaftConstants& constants,
                                                                  double lambdar) const {
  // Catches case where there is no rain present.
  if (lambdar == 0.) {
    return {0., 0.};
  }
  double dd = 2.; // micron
  std::size_t low_k = 1;
  std::size_t high_k = 10000;
  double amg_fac = 1000. * constants.pi * constants.rhow / 6.;
  double v0_numer = 0., v0_denom = 0.;
  double v3_numer = 0., v3_denom = 0.;
  for (std::size_t kk = low_k; kk != high_k + 1; ++kk) {
    double real_kk = kk;
    double dia_micron = (real_kk - 0.5) * dd;
    double dia = dia_micron * 1.e-6;
    double amg = amg_fac * pow(dia, 3);
    double vt = 0;
    if (dia_micron <= 134.43) {
      vt = 4.5795e3 * pow(amg, 2./3.);
    } else if (dia_micron <= 1511.64) {
      vt = 4.962e1 * pow(amg, 1./3.);
    } else if (dia_micron <= 3477.84) {
      vt = 1.732e1 * pow(amg, 1./6.);
    } else {
      vt = 9.17;
    }
    v0_numer += vt * exp(-lambdar * dia);
    v0_denom += exp(-lambdar * dia);
    v3_numer += vt * pow(dia, 3.) * exp(-lambdar * dia);
    v3_denom += pow(dia, 3.) * exp(-lambdar * dia);
  }
  v0_numer *= dd * 1.e-6;
  v0_denom *= dd * 1.e-6;
  v3_numer *= dd * 1.e-6;
  v3_denom *= dd * 1.e-6;
  if (v0_denom < 1.e-30) {
    v0_denom = 1.e-30;
  }
  if (v3_denom < 1.e-30) {
    v3_denom = 1.e-30;
  }
  // Number and mass fall speeds in m/s.
  double v0 = v0_numer / v0_denom;
  double v3 = v3_numer / v3_denom;
  return {v0, v3};
}

double Sedimentation::rho_fac(const RainshaftConstants& constants, double rho_dry) const {
  return pow(1.e5 / (rho_dry * constants.rdry * 273.15), 0.54);
}
