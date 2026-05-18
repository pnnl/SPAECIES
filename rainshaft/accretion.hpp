#ifndef ACCRETION_HPP
#define ACCRETION_HPP
#include <vector>
#include <cmath>

#include "rainshaft_process.hpp"

class Accretion : public RainshaftProcess {

private:
  const double prefactor;
  const double qc_exponent;
  const double qr_exponent;

public:

  Accretion(const double prefactor, const double qc_exponent,
            const double qr_exponent);

  // Calculate tendency from current state.
  void calc_tend(const RainshaftConstants& constants,
                 const RainshaftGrid& grid,
                 const StateConst& state,
                 const RainshaftDerivedVars& dvars,
                 Tendency& tend) const;

  void calc_tend_jac(const RainshaftConstants &constants,
                             const RainshaftGrid &grid,
                             const StateConst& state,
                             const RainshaftDerivedVars &dvars,
                             Matrix jac) const;

};

#endif // ACCRETION_HPP
