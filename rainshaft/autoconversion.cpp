#include "autoconversion.hpp"
#include "derivatives.hpp"
#include <stdexcept>

double Autoconversion::drop_radius_to_inv_mass(const RainshaftConstants& constants,
                                               const double new_rain_radius) {
  return 3. / (4. * constants.pi * constants.rhow
               * new_rain_radius * new_rain_radius * new_rain_radius);
}

Autoconversion::Autoconversion(const RainshaftConstants& constants,
                               const double prefactor, const double nc_exponent,
                               const double qc_exponent, const double new_rain_radius)
  : prefactor(prefactor), nc_exponent(nc_exponent), qc_exponent(qc_exponent),
    new_rain_inv_mass(drop_radius_to_inv_mass(constants, new_rain_radius))
{
}

void Autoconversion::calc_tend(const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const StateConst& state,
                              const RainshaftDerivedVars& dvars,
                              Tendency& tend) const {
  VarConst nc = state.get_variable("nc");
  VarConst qc = state.get_variable("qc");
  VarMut nc_tend = tend.get_variable("nc_tend");
  VarMut qc_tend = tend.get_variable("qc_tend");
  VarMut nr_tend = tend.get_variable("nr_tend");
  VarMut qr_tend = tend.get_variable("qr_tend");
  double qr_auto_tend;
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    if (qc[il] >= 1.e-8) {
      qr_auto_tend = prefactor * pow(qc[il], qc_exponent)
        * pow(nc[il] * 1.e-6 * dvars.rho_dry[il], nc_exponent);
      nc_tend[il] -= qr_auto_tend * nc[il] / qc[il];
      qc_tend[il] -= qr_auto_tend;
      nr_tend[il] += qr_auto_tend * new_rain_inv_mass;
      qr_tend[il] += qr_auto_tend;
    }
  }
}

void Autoconversion::calc_tend_jac(const RainshaftConstants &constants,
                                  const RainshaftGrid &grid,
                                  const StateConst& state,
                                  const RainshaftDerivedVars &dvars,
                                  Matrix jac) const
{
  // Not implemented
}
