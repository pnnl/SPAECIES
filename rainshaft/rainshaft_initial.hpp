#ifndef RAINSHAFT_INITIAL_HPP
#define RAINSHAFT_INITIAL_HPP
#include "rainshaft_constants.hpp"
#include "rainshaft_derived_vars.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_types.hpp"
#include "saturation.hpp"

bool warm_adiabatic_initial_condition(const RainshaftConstants &constants,
                                      const RainshaftGrid &grid,
                                      const SaturationFormulae &sat_form,
                                      double srf_temp,
                                      double lapse_rate,
                                      double rel_hum_init,
                                      State& initial_state);

#endif // RAINSHAFT_INITIAL_HPP
