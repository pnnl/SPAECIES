#ifndef EXPLICIT_INTEGRATOR_HPP
#define EXPLICIT_INTEGRATOR_HPP
#include "sundials_integrator.hpp"

class ExplicitIntegrator : public SundialsIntegrator {

public:

  double dt;
  
  ExplicitIntegrator(const RainshaftConstants* constants,
                     const double dt_in,
                     const RainshaftGrid* grid,
                     const RainshaftProcess* process_exp,
                     const RainshaftProcess* process_imp,
                     sundials::Context *sun_ctxt);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const RainshaftState& initial_state) const;

};

#endif // EXPLICIT_INTEGRATOR_HPP
