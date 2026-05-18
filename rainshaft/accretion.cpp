#include "accretion.hpp"
#include "derivatives.hpp"
#include <stdexcept>

Accretion::Accretion(const double prefactor, const double qc_exponent,
                               const double qr_exponent)
  : prefactor(prefactor), qc_exponent(qc_exponent), qr_exponent(qr_exponent)
{
}

void Accretion::calc_tend(const RainshaftConstants& constants,
                          const RainshaftGrid& grid,
                          const StateConst& state,
                          const RainshaftDerivedVars& dvars,
                          Tendency& tend) const {
  VarConst nc = state.get_variable("nc");
  VarConst qc = state.get_variable("qc");
  VarConst qr = state.get_variable("qr");
  VarMut nc_tend = tend.get_variable("nc_tend");
  VarMut qc_tend = tend.get_variable("qc_tend");
  VarMut qr_tend = tend.get_variable("qr_tend");
  double qr_accr_tend;
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    if (qc[il] >= constants.qsmall && qr[il] >= constants.qsmall) {
      qr_accr_tend = prefactor * pow(qc[il], qc_exponent)
        * pow(qr[il], qr_exponent);
      nc_tend[il] -= qr_accr_tend * nc[il] / qc[il];
      qc_tend[il] -= qr_accr_tend;
      qr_tend[il] += qr_accr_tend;
    }
  }
}

void Accretion::calc_tend_jac(const RainshaftConstants &constants,
                              const RainshaftGrid &grid,
                              const StateConst& state,
                              const RainshaftDerivedVars &dvars,
                              Matrix jac) const
{
  // Not implemented
}
