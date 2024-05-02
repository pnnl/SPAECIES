#ifndef EXPLICIT_INTEGRATOR_HPP
#define EXPLICIT_INTEGRATOR_HPP
#include "sundials_integrator.hpp"

class ExplicitIntegrator : public SundialsIntegrator<1> {

public:

  ExplicitIntegrator(const RainshaftConstants& constants,
                     const RainshaftGrid& grid,
                     const RainshaftProcess* const process,
                     const int steps_per_output = -1);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const RainshaftState& initial_state) const;

};

#endif // EXPLICIT_INTEGRATOR_HPP
