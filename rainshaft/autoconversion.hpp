#ifndef AUTOCONVERSION_HPP
#define AUTOCONVERSION_HPP
#include <vector>
#include <cmath>

#include "rainshaft_process.hpp"

class Autoconversion : public RainshaftProcess {

private:
  const double prefactor;
  const double nc_exponent;
  const double qc_exponent;
  const double new_rain_inv_mass; // inverse mass of new raindrops, in 1/kg

  double drop_radius_to_inv_mass(const RainshaftConstants& constants, const double radius);

public:

   // new_rain_radius is the average size of newly formed raindrops, in meters
  Autoconversion(const RainshaftConstants& constants,
                 const double prefactor, const double nc_exponent,
                 const double qc_exponent, const double new_rain_radius);

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

#endif // AUTOCONVERSION_HPP
