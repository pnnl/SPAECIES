#ifndef SEDIMENTATION_HPP
#define SEDIMENTATION_HPP
#include <vector>
#include <optional>
#include <tuple>

#include "spaecies.hpp"

#include "lookup_linear.hpp"
#include "rainshaft_process.hpp"

class Sedimentation : public RainshaftProcess {

public:
  using Speeds = std::tuple<double, double>;

  // Constructor allowing use of lookup table.
  Sedimentation(const RainshaftConstants& constants, bool use_v_table,
                bool use_numerical_integration);

  // Calculate tendency from current state.
  RainshaftTendency calc_tend(const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const spaecies::VariableArrayView<double>& state,
                              const RainshaftDerivedVars& dvars) const;

  // For a given value of lambdar, what are the rain number and mass fall speeds?
  Speeds rain_fall_speeds(const RainshaftConstants& constants,
                                       double rho, double lambdar) const;

  // What would the speeds for a given lambdar be at standard temperature and
  // pressure?
  Speeds rain_fall_speeds_stp(const RainshaftConstants& constants,
                                           double lambdar) const;

  // What would the speeds for a given lambdar be at standard temperature and
  // pressure?
  // This version ignores any lookup table, if present, and always returns fall
  // speeds calculated with incomplete gamma functions.
  Speeds rain_fall_speeds_stp_gamma(const RainshaftConstants& constants,
                                                 double lambdar) const;

  // What would the speeds for a given lambdar be at standard temperature and
  // pressure?
  // This version ignores any lookup table, if present, and always returns fall
  // speeds calculated with numerical integration.
  Speeds rain_fall_speeds_stp_numerical(const RainshaftConstants& constants,
                                                     double lambdar) const;

  // Conversion factor from STP fall speed to actual fall speed given air density.
  double rho_fac(const RainshaftConstants& constants,
                 double rho_dry) const;

  const bool use_numerical_integration;

private:

  // Particle velocity lookup tables.
  // This is currently hard-coded to P3 settings, i.e. it contains 20 entries
  // for values every 10 microns between 5 and 195 micron, and 280 entries for
  // values every 30 microns between 195 and 8595 micron.
  std::optional<LookupLinear> v0_table;
  std::optional<LookupLinear> v3_table;

};

#endif // SEDIMENTATION_HPP
