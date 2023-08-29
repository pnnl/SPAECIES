#ifndef EVAPORATION_HPP
#define EVAPORATION_HPP
#include <vector>

#include "rainshaft_process.hpp"
#include "saturation.hpp"

class Evaporation : public RainshaftProcess {

public:

  // SPS: Should use a different pointer type to guarantee this
  // does not outlast the sat_form object.
  Evaporation(const SaturationFormulae* sat_form_in);

  // Calculate tendency from current state.
  RainshaftTendency calc_tend(const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const RainshaftState& state,
                              const RainshaftDerivedVars& dvars) const;

private:

  const SaturationFormulae *sat_form;

};

#endif // EVAPORATION_HPP
