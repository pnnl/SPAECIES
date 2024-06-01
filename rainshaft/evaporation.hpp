#ifndef EVAPORATION_HPP
#define EVAPORATION_HPP
#include <vector>
#include <optional>

#include "lookup_linear.hpp"
#include "rainshaft_process.hpp"
#include "saturation.hpp"

class Evaporation : public RainshaftProcess
{

public:
  // SPS: Should use a different pointer type to guarantee this
  // does not outlast the sat_form object.
  Evaporation(const RainshaftConstants &constants,
              const SaturationFormulae *sat_form_in, bool use_v_table,
              bool use_numerical_integration);

  // Calculate tendency from current state.
  RainshaftTendency calc_tend(const RainshaftConstants &constants,
                              const RainshaftGrid &grid,
                              const RainshaftState &state,
                              const RainshaftDerivedVars &dvars) const;

  void calc_tend_jac_prod(const RainshaftConstants &constants,
                          const RainshaftGrid &grid,
                          const RainshaftState &state,
                          const RainshaftDerivedVars &dvars,
                          const double *const vec,
                          double *const prod) const;

  void calc_tend_jac(const RainshaftConstants &constants,
                             const RainshaftGrid &grid,
                             const RainshaftState &state,
                             const RainshaftDerivedVars &dvars,
                             SUNMatrix jac) const;

  // Calculate characteristic velocity used for velocity calculation.
  double calc_v_evap(const RainshaftConstants &constants, double lambdar) const;

  // Calculate characteristic velocity used for velocity calculation.
  // This version ignores the lookup table, if present, and always just
  // calculates using incomplete gamma functions.
  double calc_v_evap_gamma(const RainshaftConstants &constants, double lambdar) const;

  // Calculate characteristic velocity using numerical integration.
  double calc_v_evap_numerical(const RainshaftConstants &constants, double lambdar) const;

  const bool use_numerical_integration;

private:
  // Water vapor saturation formulae.
  const SaturationFormulae *sat_form;

  // Particle velocity lookup table.
  // This is currently hard-coded to P3 settings, i.e. it contains 20 entries
  // for values every 10 microns between 5 and 195 micron, and 280 entries for
  // values every 30 microns between 195 and 8595 micron.
  std::optional<LookupLinear> v_table;
};

#endif // EVAPORATION_HPP
