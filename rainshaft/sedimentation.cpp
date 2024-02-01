#include "sedimentation.hpp"
#include <cstddef>
#include <cmath>
#include <iostream>
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
      auto speeds = use_numerical_integration ?
        rain_fall_speeds_stp_numerical(constants, lambdar)
        :  rain_fall_speeds_stp_gamma(constants, lambdar);
      v0_values[i] = speeds[0];
      v3_values[i] = speeds[1];
    }
    v0_table.emplace(range_bounds, spacings, v0_values);
    v3_table.emplace(range_bounds, spacings, v3_values);
  }
}

RainshaftTendency Sedimentation::calc_tend(const RainshaftConstants& constants,
                                           const RainshaftGrid& grid,
                                           const RainshaftState& state,
                                           const RainshaftDerivedVars& dvars) const {
  std::vector<double> t_tend, q_tend, nr_tend, qr_tend;
  double nr_tend_lev, qr_tend_lev;
  std::vector<double> v0, v3;
  t_tend.push_back(0.);
  q_tend.push_back(0.);
  // SPS: Should have a utility rather than duplicating lambdar calculation.
  double lambdar_top = constants.pi * constants.rhow * constants.nr_top / constants.qr_top;
  lambdar_top = cbrt(lambdar_top);
  std::vector<double> speeds_top = rain_fall_speeds(constants, constants.rho_top, lambdar_top);
  double v0_top = speeds_top[0], v3_top = speeds_top[3];
  std::vector<double> speeds = rain_fall_speeds(constants, dvars.rho_dry[0], dvars.lambdar[0]);
  v0.push_back(speeds[0]);
  v3.push_back(speeds[3]);
  nr_tend_lev = v0_top * constants.nr_top * constants.rho_top - v0[0]*state.nr[0]*dvars.rho_dry[0];
  nr_tend_lev /= dvars.dz[0] * dvars.rho_dry[0];
  qr_tend_lev = v3_top * constants.qr_top * constants.rho_top - v3[0]*state.qr[0]*dvars.rho_dry[0];
  qr_tend_lev /= dvars.dz[0] * dvars.rho_dry[0];
  nr_tend.push_back(nr_tend_lev);
  qr_tend.push_back(qr_tend_lev);
  for (std::size_t il = 1; il != grid.nlev; ++il) {
    t_tend.push_back(0.);
    q_tend.push_back(0.);
    speeds = rain_fall_speeds(constants, dvars.rho_dry[il], dvars.lambdar[il]);
    v0.push_back(speeds[0]);
    v3.push_back(speeds[3]);
    nr_tend_lev = v0[il-1]*state.nr[il-1]*dvars.rho_dry[il-1] - v0[il]*state.nr[il]*dvars.rho_dry[il];
    nr_tend_lev /= dvars.dz[il] * dvars.rho_dry[il];
    qr_tend_lev = v3[il-1]*state.qr[il-1]*dvars.rho_dry[il-1] - v3[il]*state.qr[il]*dvars.rho_dry[il];
    qr_tend_lev /= dvars.dz[il] * dvars.rho_dry[il];

    nr_tend.push_back(nr_tend_lev);
    qr_tend.push_back(qr_tend_lev);
  }
  return RainshaftTendency(t_tend, q_tend, nr_tend, qr_tend);
}


// Jacobian of sedimentation
RainshaftTendencyJac Sedimentation::calc_tend_jac(const RainshaftConstants& constants,
                                           const RainshaftGrid& grid,
                                           const RainshaftState& state,
                                           const RainshaftDerivedVars& dvars) const {
  double* t_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* q_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* nr_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* qr_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};

  std::vector<double> v0, v1, v2, v3, v4;
  std::vector<double> speeds;

  speeds = rain_fall_speeds(constants, dvars.rho_dry[0], dvars.lambdar[0]);
  // std::cout << speeds[0] << " " << speeds[1] << " " << speeds[2] << " " << speeds[3] << " " << speeds[4] << std::endl;


  // d(f_nr)/d(nr)
  for (std::size_t il = 2*grid.nlev; il < 3*grid.nlev; ++il) {
    for (std::size_t jl = 2*grid.nlev; jl < 3*grid.nlev; ++jl) {
      if (il == jl) {
        // std::cout << dvars.lambdar[jl-2*grid.nlev] << std::endl;
        speeds = rain_fall_speeds(constants, dvars.rho_dry[jl-2*grid.nlev], dvars.lambdar[jl-2*grid.nlev]);
        nr_tend_jac[jl*4*grid.nlev + il] = -speeds[0]*dvars.rho_dry[jl-2*grid.nlev] - (speeds[0] - speeds[1])/3.0*dvars.rho_dry[jl-2*grid.nlev];
        nr_tend_jac[jl*4*grid.nlev + il] /= dvars.dz[jl-2*grid.nlev] * dvars.rho_dry[jl-2*grid.nlev];
        // std::cout << -speeds[0]*dvars.rho_dry[jl-2*grid.nlev] - (speeds[0] - speeds[1])/3.0*dvars.rho_dry[jl-2*grid.nlev] << std::endl;
        // std::cout << speeds[0] << std::endl;
      }

      if (jl == il-1) {
        speeds = rain_fall_speeds(constants, dvars.rho_dry[jl-2*grid.nlev], dvars.lambdar[jl-2*grid.nlev]);
        nr_tend_jac[jl*4*grid.nlev + il] = speeds[0]*dvars.rho_dry[jl-2*grid.nlev] + (speeds[0] - speeds[1])/3.0*dvars.rho_dry[jl-2*grid.nlev];
        nr_tend_jac[jl*4*grid.nlev + il] /= dvars.dz[jl-2*grid.nlev] * dvars.rho_dry[jl-2*grid.nlev];
      }
    }
  }

  // d(f_nr)/d(qr)
  for (std::size_t il = 2*grid.nlev; il < 3*grid.nlev; ++il) {
    for (std::size_t jl = 3*grid.nlev; jl < 4*grid.nlev; ++jl) {
      if (il+grid.nlev == jl) {
        speeds = rain_fall_speeds(constants, dvars.rho_dry[jl-3*grid.nlev], dvars.lambdar[jl-3*grid.nlev]);
        nr_tend_jac[jl*4*grid.nlev + il] = -(speeds[0] - speeds[1])*-pow(dvars.lambdar[jl-3*grid.nlev], 3)/(3.0*M_PI*constants.rhow)*dvars.rho_dry[jl-3*grid.nlev];
        nr_tend_jac[jl*4*grid.nlev + il] /= dvars.dz[jl-3*grid.nlev] * dvars.rho_dry[jl-3*grid.nlev];
        // std::cout << -(speeds[0] - speeds[1])*pow(dvars.lambdar[jl-3*grid.nlev], 3)/3.0*M_PI*constants.rhow*dvars.rho_dry[jl-3*grid.nlev] << std::endl;
      }

      if (jl == il-1 + grid.nlev) {
        speeds = rain_fall_speeds(constants, dvars.rho_dry[jl-3*grid.nlev], dvars.lambdar[jl-3*grid.nlev]);
        nr_tend_jac[jl*4*grid.nlev + il] = (speeds[0] - speeds[1])*-pow(dvars.lambdar[jl-3*grid.nlev], 3)/(3.0*M_PI*constants.rhow)*dvars.rho_dry[jl-3*grid.nlev];
        nr_tend_jac[jl*4*grid.nlev + il] /= dvars.dz[jl-3*grid.nlev] * dvars.rho_dry[jl-3*grid.nlev];
      }
    }
  }

  // d(f_qr)/d(nr)
  for (std::size_t il = 3*grid.nlev; il < 4*grid.nlev; ++il) {
    for (std::size_t jl = 2*grid.nlev; jl < 3*grid.nlev; ++jl) {
      if (il == jl+grid.nlev) {
        // std::cout << "switching " << speeds[3] << " " << speeds[4] << std::endl;
        
        if (dvars.lambdar[jl-2*grid.nlev] == 0.0) {
          speeds = rain_fall_speeds(constants, dvars.rho_dry[jl-2*grid.nlev], 1.e5);
          qr_tend_jac[jl*4*grid.nlev + il] = -(speeds[3] - speeds[4])*4.0*(M_PI*constants.rhow/(3.0*pow(1.e5, 3)))*dvars.rho_dry[jl-2*grid.nlev];
        } else {
          speeds = rain_fall_speeds(constants, dvars.rho_dry[jl-2*grid.nlev], dvars.lambdar[jl-2*grid.nlev]);
          qr_tend_jac[jl*4*grid.nlev + il] = -(speeds[3] - speeds[4])*4.0*(M_PI*constants.rhow/(3.0*pow(dvars.lambdar[jl-2*grid.nlev], 3)))*dvars.rho_dry[jl-2*grid.nlev];
        }
        qr_tend_jac[jl*4*grid.nlev + il] /= dvars.dz[jl-2*grid.nlev] * dvars.rho_dry[jl-2*grid.nlev];
      }

      if (jl+grid.nlev == il-1) {
        // speeds = rain_fall_speeds(constants, dvars.rho_dry[jl-2*grid.nlev], dvars.lambdar[jl-2*grid.nlev]);
        if (dvars.lambdar[jl-2*grid.nlev] == 0.0) {
          speeds = rain_fall_speeds(constants, dvars.rho_dry[jl-2*grid.nlev], 1.e5);
          qr_tend_jac[jl*4*grid.nlev + il] = (speeds[3] - speeds[4])*4.0*(M_PI*constants.rhow/(3.0*pow(1.e5, 3)))*dvars.rho_dry[jl-2*grid.nlev];
          
        } else {
          speeds = rain_fall_speeds(constants, dvars.rho_dry[jl-2*grid.nlev], dvars.lambdar[jl-2*grid.nlev]);
          qr_tend_jac[jl*4*grid.nlev + il] = (speeds[3] - speeds[4])*4.0*(M_PI*constants.rhow/(3.0*pow(dvars.lambdar[jl-2*grid.nlev], 3)))*dvars.rho_dry[jl-2*grid.nlev];
        }
        qr_tend_jac[jl*4*grid.nlev + il] /= dvars.dz[jl-2*grid.nlev] * dvars.rho_dry[jl-2*grid.nlev];
      }
    }
  }

  // d(f_qr)/d(qr)
  // v3[il-1]*state.qr[il-1]*dvars.rho_dry[il-1] - v3[il]*state.qr[il]*dvars.rho_dry[il];
  for (std::size_t il = 3*grid.nlev; il < 4*grid.nlev; ++il) {
    for (std::size_t jl = 3*grid.nlev; jl < 4*grid.nlev; ++jl) {
      if (il == jl) {
        speeds = rain_fall_speeds(constants, dvars.rho_dry[jl-3*grid.nlev], dvars.lambdar[jl-3*grid.nlev]);
        qr_tend_jac[jl*4*grid.nlev + il] = -speeds[3]*dvars.rho_dry[jl-3*grid.nlev] - (speeds[3] - speeds[4])*-4.0/3.0*dvars.rho_dry[jl-3*grid.nlev];
        qr_tend_jac[jl*4*grid.nlev + il] /= dvars.dz[jl-3*grid.nlev] * dvars.rho_dry[jl-3*grid.nlev];
      }

      if (jl == il-1) {
        speeds = rain_fall_speeds(constants, dvars.rho_dry[jl-3*grid.nlev], dvars.lambdar[jl-3*grid.nlev]);
        qr_tend_jac[jl*4*grid.nlev + il] = speeds[3]*dvars.rho_dry[jl-3*grid.nlev] + (speeds[3] - speeds[4])*-4.0/3.0*dvars.rho_dry[jl-3*grid.nlev];
        qr_tend_jac[jl*4*grid.nlev + il] /= dvars.dz[jl-3*grid.nlev] * dvars.rho_dry[jl-3*grid.nlev];
      }
    }
  }


  // int nz = 4*grid.nlev;
  // for (int i = nz/2; i != nz; ++i) {
  //   for (int j = nz/2; j != nz; ++j) {
  //     std::cout << nr_tend_jac[i*nz + j] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // std::cout << std::endl;




  return RainshaftTendencyJac(t_tend_jac, q_tend_jac, nr_tend_jac, qr_tend_jac);
}


std::vector<double> Sedimentation::rain_fall_speeds(const RainshaftConstants& constants,
                                                    double rho_dry, double lambdar) const {
  // Catches case where there is no rain present.
  if (lambdar == 0.) {
    std::vector<double> speeds = {0., 0., 0., 0., 0.};
    return speeds;
  }
  auto speeds = rain_fall_speeds_stp(constants, lambdar);
  // Calculate correction for air density, with reference state at 1000 hPa and 273.15 K.
  double rf = rho_fac(constants, rho_dry);
  speeds[0] *= rf;
  speeds[1] *= rf;
  speeds[2] *= rf;
  speeds[3] *= rf;
  speeds[4] *= rf;
  return speeds;
}

std::vector<double> Sedimentation::rain_fall_speeds_stp(const RainshaftConstants& constants,
                                                        double lambdar) const {
  if (v0_table.has_value()) {
    double d_micron = 1.e6 / lambdar;
    double v0 = v0_table->lookup_value(d_micron);
    double v3 = v3_table->lookup_value(d_micron);
    return std::vector<double>{v0, v3};
  } else {
    if (use_numerical_integration) {
      return rain_fall_speeds_stp_numerical(constants, lambdar);
    } else {
      return rain_fall_speeds_stp_gamma(constants, lambdar);
    }
  }
}

std::vector<double> Sedimentation::rain_fall_speeds_stp_gamma(const RainshaftConstants& constants,
                                                              double lambdar) const {
  // Catches case where there is no rain present.
  if (lambdar == 0.) {
    std::vector<double> speeds = {0., 0., 0., 0., 0.};
    return speeds;
  }
  // Number and mass fall speeds in m/s.
  double v0(0.), v1(0.), v2(0.), v3(0.), v4(0.);
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


  // Start with number term for D <= 134.43 micron.
  v1 += 4579.5 * d2g_2third * tgamma_lower(4, lambdar * 1.3443e-4) / (lambdar * lambdar);
  // Number term for 134.43 micron < D <= 1511.64 micron.
  v1 += 49.62 * d2g_1third
    * (tgamma(3, lambdar * 1.3443e-4) - tgamma(3, lambdar * 1.51164e-3)) / lambdar;
  // Number term for 1511.64 micron < D <= 3477.84 micron.
  v1 += 17.32 * d2g_1sixth
    * (tgamma(2.5, lambdar * 1.51164e-3) - tgamma(2.5, lambdar * 3.47784e-3)) / sqrt(lambdar);
  // Number term for 3477.84 micron < D.
  v1 += 9.17 * tgamma(2, lambdar * 3.47784e-3);


  // Start with number term for D <= 134.43 micron.
  v2 += 4579.5 * d2g_2third * tgamma_lower(5, lambdar * 1.3443e-4) / (lambdar * lambdar);
  // Number term for 134.43 micron < D <= 1511.64 micron.
  v2 += 49.62 * d2g_1third
    * (tgamma(4, lambdar * 1.3443e-4) - tgamma(4, lambdar * 1.51164e-3)) / lambdar;
  // Number term for 1511.64 micron < D <= 3477.84 micron.
  v2 += 17.32 * d2g_1sixth
    * (tgamma(3.5, lambdar * 1.51164e-3) - tgamma(3.5, lambdar * 3.47784e-3)) / sqrt(lambdar);
  // Number term for 3477.84 micron < D.
  v2 += 9.17 * tgamma(3, lambdar * 3.47784e-3);
  // Include the normalization gamma(2) for v2.
  v2 /= 2.;


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


  // Mass term for D <= 134.43 micron.
  v4 += 4579.5 * d2g_2third * tgamma_lower(7, lambdar * 1.3443e-4) / (lambdar * lambdar);
  // Mass term for 134.43 micron < D <= 1511.64 micron.
  v4 += 49.62 * d2g_1third
    * (tgamma(6, lambdar * 1.3443e-4) - tgamma(6, lambdar * 1.51164e-3)) / lambdar;
  // Mass term for 1511.64 micron < D <= 3477.84 micron.
  v4 += 17.32 * d2g_1sixth
    * (tgamma(5.5, lambdar * 1.51164e-3) - tgamma(5.5, lambdar * 3.47784e-3)) / sqrt(lambdar);
  // Mass term for 3477.84 micron < D.
  v4 += 9.17 * tgamma(5, lambdar * 3.47784e-3);
  // Include the normalization gamma(3) for v3.
  v4 /= 24.;

  std::vector<double> speeds = {v0, v1, v2, v3, v4};
  return speeds;
}

std::vector<double> Sedimentation::rain_fall_speeds_stp_numerical(const RainshaftConstants& constants,
                                                                  double lambdar) const {
  // Catches case where there is no rain present.
  if (lambdar == 0.) {
    std::vector<double> speeds = {0., 0.};
    return speeds;
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
  std::vector<double> speeds = {v0, v3};
  return speeds;
}

double Sedimentation::rho_fac(const RainshaftConstants& constants, double rho_dry) const {
  return pow(1.e5 / (rho_dry * constants.rdry * 273.15), 0.54);
}
