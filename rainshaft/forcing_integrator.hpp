#ifndef FORCING_INTEGRATOR_HPP
#define FORCING_INTEGRATOR_HPP

#include "sundials_integrator.hpp"

class ForcingIntegrator : public SundialsIntegrator<2>
{
private:
  const double dt;
  const double dt_unforced;
  const double dt_forced;
  const bool postprocess;

public:
  ForcingIntegrator(const RainshaftConstants &constants,
                    const RainshaftGrid &grid,
                    const RainshaftProcess *const process_forced,
                    const RainshaftProcess *const process_unforced,
                    const VarDescList &state_descs,
                    const VarDescList &tend_descs,
                    const double dt,
                    const double dt_unforced = 0,
                    const double dt_forced = 0,
                    const bool postprocess = false,
                    const int steps_per_output = -1);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const StateConst &initial_state) const;
};

#endif // FORCING_INTEGRATOR_HPP
