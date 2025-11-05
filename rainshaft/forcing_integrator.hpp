#ifndef FORCING_INTEGRATOR_HPP
#define FORCING_INTEGRATOR_HPP

#include "sundials_integrator.hpp"
#include "sedimentation.hpp"

class ForcingIntegrator : public SundialsIntegrator<2>
{
private:
  const double dt;
  const double dt_forced;
  const double dt_unforced;
  const bool cfl_substep;
  const bool postprocess;

public:
  ForcingIntegrator(const RainshaftConstants &constants,
                    const RainshaftGrid &grid,
                    const Sedimentation *const process_forced,
                    const RainshaftProcess *const process_unforced,
                    const VarDescList &state_descs,
                    const VarDescList &tend_descs,
                    const double dt,
                    const double dt_forced = 0,
                    const double dt_unforced = 0,
                    const bool cfl_substep = false,
                    const bool postprocess = false,
                    const bool regularize_lambdar = true,
                    const int steps_per_output = -1);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const StateConst &initial_state,
                              int& error_flag) const;
};

#endif // FORCING_INTEGRATOR_HPP
