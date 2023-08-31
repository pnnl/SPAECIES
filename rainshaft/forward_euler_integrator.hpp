#ifndef FORWARD_EULER_INTEGRATOR_HPP
#define FORWARD_EULER_INTEGRATOR_HPP
#include "rainshaft_integrator.hpp"

class ForwardEulerIntegrator : public RainshaftIntegrator {

public:

  // Constructor requires timestep.
  ForwardEulerIntegrator(const RainshaftConstants* constants,
                         const RainshaftGrid* grid,
                         const RainshaftProcess* process,
                         sundials::Context *sun_ctxt, double dt_in);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const RainshaftState& initial_state);

  const double dt;
};

#endif // FORWARD_EULER_INTEGRATOR_HPP
