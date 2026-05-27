#include "cloud_sedimentation.hpp"
#include <cstddef>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
using boost::math::tgamma, boost::math::tgamma_lower;
using std::pow, std::sqrt, std::cbrt, std::exp;

CloudSedimentation::CloudSedimentation()
{}

double CloudSedimentation::calc_moment_flux(const double moment,
                                            const double rho_dry,
                                            const double v) const
{
  return v * moment * rho_dry;
}

double CloudSedimentation::calc_moment_tend(const double dz,
                                            const double rho_dry,
                                            const double flux_bot,
                                            const double flux_top) const
{
  const double air_mass_per_area = dz * rho_dry;
  return (flux_top - flux_bot) / air_mass_per_area;
}

double CloudSedimentation::calc_muc(double nc, double rho_dry) const
{
  // Convert nc to #/cm^3.
  double nc_cgs = std::max(nc, 1.e-16) * 1.e-6 * rho_dry; // NSMALL=10^-16 from EAMxx
  double denom_fac = 0.0005714 * nc_cgs + 0.2714;
  double muc = (1. / (denom_fac*denom_fac)) - 1.;
  // Output limited to the range [2, 15].
  return std::min(std::max(2., muc), 15.);
}

double CloudSedimentation::calc_lambdac(const RainshaftConstants& constants,
                                        double nc, double qc, double muc) const
{
  if (qc == 0.) {
    // Default value when no mass present.
    return 0.;
  }
  double muc_poly = constants.pi * constants.rhow * (muc + 3.) * (muc + 2.) * (muc + 1.) / 6.;
  return std::cbrt(muc_poly * nc / qc);
}

// For a given size distribution+temperature, what are the cloud number and mass fall speeds?
std::tuple<double, double> CloudSedimentation::cloud_fall_speeds(const RainshaftConstants &constants, 
                                                                 double t, double lambdac, double muc) const
{
  if (lambdac == 0.) {
    return {0., 0.};
  }
  // Approximate dynamic viscosity based on temperature of the air.
  const double viscosity = 1.496e-6 * std::pow(t, 1.5) / (t + 120.);
  // Constant (w.r.t. particle diameter) factor used in Stokes' law.
  const double acn = constants.g * constants.rhow / (18. * viscosity);

  // This factor theoretically could change due to changes in parameterization assumptions/tuning.
  // But in practice, there are some reasons it should be about 2 for gamma-distributed cloud
  // and it may be worth considering optimizing for integer values in the future to avoid
  // the pow/gamma calls.
  const double bcn = 2.;

  // Fall speeds for 0th and 3rd moments (number and mass) based on integrating Stokes' law
  // analytically over gamma-distributed diameter.
  const double lam_fac = acn * std::pow(lambdac, -bcn);
  const double v0 = lam_fac * tgamma(muc + 1 + bcn) / tgamma(muc + 1);
  const double v3 = lam_fac * tgamma(muc + 4 + bcn) / tgamma(muc + 4);

  return {v0, v3};
}

void CloudSedimentation::calc_tend(const RainshaftConstants &constants,
                                   const RainshaftGrid &grid,
                                   const StateConst &state,
                                   const RainshaftDerivedVars &dvars,
                                   Tendency& tend) const
{
  VarConst t = state.get_variable("T");
  VarConst nc = state.get_variable("nc");
  VarConst qc = state.get_variable("qc");
  VarMut nc_tend = tend.get_variable("nc_tend");
  VarMut qc_tend = tend.get_variable("qc_tend");

  const double muc_top = calc_muc(constants.nc_top, constants.rho_top);
  const double lambdac_top = calc_lambdac(constants, constants.nc_top, constants.qc_top, muc_top);

  // Note that ideally we should probably have a t_top boundary condition here, rather than use
  // the top level temperature.
  auto [v0_prev, v3_prev] = cloud_fall_speeds(constants, t[0], lambdac_top, muc_top);
  double nc_prev = constants.nc_top;
  double qc_prev = constants.qc_top;
  double nc_flux_prev = calc_moment_flux(nc_prev, constants.rho_top, v0_prev);
  double qc_flux_prev = calc_moment_flux(qc_prev, constants.rho_top, v3_prev);

  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const double muc = calc_muc(nc[il], dvars.rho_dry[il]);
    const double lambdac = calc_lambdac(constants, nc[il], qc[il], muc);
    const auto [v0, v3] = cloud_fall_speeds(constants, t[il], lambdac, muc);
    const double nc_flux = calc_moment_flux(nc[il], dvars.rho_dry[il], v0);
    const double qc_flux = calc_moment_flux(qc[il], dvars.rho_dry[il], v3);
    
    nc_tend[il] += calc_moment_tend(dvars.dz[il], dvars.rho_dry[il], nc_flux, nc_flux_prev);
    qc_tend[il] += calc_moment_tend(dvars.dz[il], dvars.rho_dry[il], qc_flux, qc_flux_prev);

    nc_flux_prev = nc_flux;
    qc_flux_prev = qc_flux;
  }
}

void CloudSedimentation::calc_tend_jac(const RainshaftConstants &constants,
                                       const RainshaftGrid &grid,
                                       const StateConst& state,
                                       const RainshaftDerivedVars &dvars,
                                       Matrix jac) const
{
  // Not implemented.
}
