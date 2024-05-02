#ifndef MRI_INTEGRATOR_HPP
#define MRI_INTEGRATOR_HPP
#include "sundials_integrator.hpp"

class MRIIntegrator : public SundialsIntegrator<3>
{

private:
  const double dt_slow;
  const double multirate_ratio;
  const int order;

public:
  MRIIntegrator(const RainshaftConstants &constants,
                const RainshaftGrid &grid,
                const RainshaftProcess * const process_fast,
                const RainshaftProcess * const process_slow_exp,
                const RainshaftProcess *const process_slow_imp,
                const double dt_slow,
                const int multirate_ratio,
                const int order,
                const int steps_per_output = -1);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const RainshaftState &initial_state) const;
};

#endif // MRI_INTEGRATOR_HPP
