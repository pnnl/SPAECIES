#ifndef FORWARD_EULER_INTEGRATOR_HPP
#define FORWARD_EULER_INTEGRATOR_HPP
#include "rainshaft_integrator.hpp"

class ForwardEulerIntegrator : public RainshaftIntegrator {

public:

  // Constructor requires timestep.
  ForwardEulerIntegrator(sundials::Context *sun_ctxt, double dt_in);

  RainshaftSolution integrate(const RainshaftProcess& process,
                              double initial_time,
                              double final_time,
                              const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const RainshaftState& initial_state,
                              const RainshaftDerivedVars& initial_dvars);

  const double dt;
};

#endif // FORWARD_EULER_INTEGRATOR_HPP
