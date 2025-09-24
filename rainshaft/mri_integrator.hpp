#ifndef MRI_INTEGRATOR_HPP
#define MRI_INTEGRATOR_HPP
#include "sundials_integrator.hpp"

class MRIIntegrator : public SundialsIntegrator<3>
{

private:
  const double dt_fast;
  const double dt_slow;
  const int order;
  const double rel_tol;
  const bool postprocess;

public:
  MRIIntegrator(const RainshaftConstants &constants,
                const RainshaftGrid &grid,
                const RainshaftProcess *const process_fast,
                const RainshaftProcess *const process_slow_exp,
                const RainshaftProcess *const process_slow_imp,
                const VarDescList& state_descs,
                const VarDescList& tend_descs,
                const double dt_fast,
                const double dt_slow,
                const int order = 3,
                const double rel_tol = 1.e-4,
                const bool postprocess = false,
                const int steps_per_output = -1);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const StateConst &initial_state,
                              int& error_flag) const;
};

#endif // MRI_INTEGRATOR_HPP
