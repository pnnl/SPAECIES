#ifndef RAINSHAFT_CONSTANTS_HPP
#define RAINSHAFT_CONSTANTS_HPP

struct RainshaftConstants {
  // Pi!
  double pi;
  // Dry air gas constant (J/K/kg)
  double rdry;
  // Dry air heat capacity (J/K/kg)
  double cp;
  // Water vapor gas constant (J/K/kg)
  double rvapor;
  // Density of liquid water (kg/m^3)
  double rhow;
  // Latent heat of vaporization (J/kg)
  double latvap;
  // Water vapor/dry air molecular mass ratio (unitless)
  double epsilon_h2o;
  // Minimum hydrometeor mass mixing ratio (kg/kg)
  double qsmall;
  // Acceleration due to gravity at earth's surface (m/s^2)
  double g;
  // Minimum allowed mean rain diameter (m)
  double dr_min;
  // Maximum allowed mean rain diameter (m)
  double dr_max;
};

#endif // RAINSHAFT_CONSTANTS_HPP
