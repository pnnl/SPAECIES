#ifndef FORWARD_EULER_INTEGRATOR_HPP
#define FORWARD_EULER_INTEGRATOR_HPP
#include "sundials_integrator.hpp"

class ForwardEulerIntegrator : public SundialsIntegrator {

public:

  // Constructor requires timestep.
  ForwardEulerIntegrator(const RainshaftConstants* constants,
                         const RainshaftGrid* grid,
                         const RainshaftProcess* process);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const RainshaftState& initial_state) const;

};

#endif // FORWARD_EULER_INTEGRATOR_HPP
