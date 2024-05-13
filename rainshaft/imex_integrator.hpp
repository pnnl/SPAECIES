#ifndef IMEX_INTEGRATOR_HPP
#define IMEX_INTEGRATOR_HPP
#include "sundials_integrator.hpp"

class IMEXIntegrator : public SundialsIntegrator<2>
{
private:
  const double dt;
  const int order;

public:
  IMEXIntegrator(const RainshaftConstants &constants,
                     const RainshaftGrid &grid,
                     const RainshaftProcess *const process_exp,
                     const RainshaftProcess *const process_imp,
                     const double dt = 0,
                     const int order = 4,
                     const int steps_per_output = -1);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const RainshaftState &initial_state) const;
};

#endif // IMEX_INTEGRATOR_HPP