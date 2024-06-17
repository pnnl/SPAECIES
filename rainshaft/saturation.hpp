#ifndef SATURATION_HPP
#define SATURATION_HPP
#include "rainshaft_constants.hpp"
#include <cmath>

// SPS: Should consider making this an abstract class with different formulas,
// since there are a number of parameterizations in use and they should be
// consistent across the model.

class SaturationFormulae {

public:

  // Initialize from constants.
  SaturationFormulae(const RainshaftConstants& constants);

  // Saturation vapor pressure over liquid at a given temperature.
  double svp_liquid(double temperature) const;

  // Convert water vapor pressure to specific humidity over dry air.
  double wv_pressure_to_q_dry(double pressure_wv, double pressure_dry) const;

  // Saturation specific humidity (over mass of dry air) as a function of
  // temperature and dry pressure.
  double q_sat_dry(double temperature, double pressure_dry) const;

  // derivative of esl wrt T
  double svp_liquid_dT(double temperature) const;

  // derivative of q_sl wrt T
  double q_sat_dry_dT(double temperature, double pressure_dry) const;

private:

  // Copy of ratio of water vapor molecular mass to that of dry air.
  // Stored mainly as a convenience to not require a RainshaftConstants as an
  // argument to a bunch of functions.
  double epsilon_h2o;

};

#endif // SATURATION_HPP
