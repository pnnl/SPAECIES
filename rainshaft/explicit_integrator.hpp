#ifndef EXPLICIT_INTEGRATOR_HPP
#define EXPLICIT_INTEGRATOR_HPP
#include "rainshaft_integrator.hpp"

class ExplicitIntegrator : public RainshaftIntegrator {

public:

  ExplicitIntegrator(sundials::Context *sun_ctxt);

  RainshaftSolution integrate(const RainshaftProcess& process,
                              double initial_time,
                              double final_time,
                              const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const RainshaftState& initial_state,
                              const RainshaftDerivedVars& initial_dvars);

};

#endif // EXPLICIT_INTEGRATOR_HPP
