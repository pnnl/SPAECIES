#include "cloud_size_parameters.hpp"
#include "rainshaft_derived_vars.hpp"

double calc_muc(double nc, double rho_dry)
{
  // Convert nc to #/cm^3.
  double nc_cgs = std::max(nc, 1.e-16) * 1.e-6 * rho_dry; // NSMALL=10^-16 from EAMxx
  double denom_fac = 0.0005714 * nc_cgs + 0.2714;
  double muc = (1. / (denom_fac*denom_fac)) - 1.;
  // Output limited to the range [2, 15].
  return std::min(std::max(2., muc), 15.);
}

double calc_lambdac(const RainshaftConstants& constants,
                    double nc, double qc, double muc)
{
  if (qc == 0.) {
    // Default value when no mass present.
    return 0.;
  }
  double muc_poly = constants.pi * constants.rhow * (muc + 3.) * (muc + 2.) * (muc + 1.) / 6.;
  return std::cbrt(muc_poly * nc / qc);
}

void limit_nc(const RainshaftConstants& constants,
              const RainshaftGrid& grid,
              State& state)
{
  // Allowed diameter range of 1 to 40 microns.
  constexpr double inv_min_diameter = 1. / 1.e-6;
  constexpr double inv_max_diameter = 1. / 40.e-6;
  VarConst t = *state.get_variable("T");
  VarConst q = *state.get_variable("q");
  VarConst qc = *state.get_variable("qc");
  VarMut nc = *state.get_variable("nc");
  for (std::size_t il = 0; il!= grid.nlev; ++il) {
    // Zero nc and move on if no water mass.
    if (qc[il] <= 0.) {
      nc[il] = 0.;
      continue;
    }
    const double rho_dry = rho_dry_from_ideal_gas_law(constants.rdry,
                                                      constants.epsilon_h2o,
                                                      grid.p_mid[il],
                                                      t[il],
                                                      q[il]);
    const double muc = calc_muc(nc[il], rho_dry);
    double lambdac = calc_lambdac(constants, nc[il], qc[il], muc);
    // *Maximum* diameter sets *minimum* lambda, and vice versa.
    const double lambdac_min = (muc + 1.) * inv_max_diameter;
    const double lambdac_max = (muc + 1.) * inv_min_diameter;
    bool nc_needs_update = false;
    if (lambdac < lambdac_min) {
      lambdac = lambdac_min;
      nc_needs_update = true;
    } else if (lambdac > lambdac_max) {
      lambdac = lambdac_max;
      nc_needs_update = true;
    }
    if (nc_needs_update) {
      // Note that this doesn't exactly invert the expression for lambdac
      // in terms of nc, because it assumes muc is constant. However, it
      // may be the case that the muc limiters make muc constant in most
      // cases where this code triggers. In any case, this is the P3 code's
      // original implementation.
      nc[il] = 6 * (lambdac*lambdac*lambdac) * qc[il]
        / (constants.pi * constants.rhow * (muc + 1.) * (muc + 2.) + (muc + 3.));
    }
  }
}
