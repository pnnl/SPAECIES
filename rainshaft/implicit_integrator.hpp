#ifndef IMPLICIT_INTEGRATOR_HPP
#define IMPLICIT_INTEGRATOR_HPP
#include "sundials_integrator.hpp"

class ImplicitIntegrator : public SundialsIntegrator {

public:

  double dt;

  ImplicitIntegrator(const RainshaftConstants* constants,
                     const double dt_in,
                     const RainshaftGrid* grid,
                     const RainshaftProcess* process,
                     sundials::Context *sun_ctxt);



  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const RainshaftState& initial_state) const;

};

#endif // IMPLICIT_INTEGRATOR_HPP
