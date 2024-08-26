#ifndef FIXED_SUBSTEP_INTEGRATOR_HPP
#define FIXED_SUBSTEP_INTEGRATOR_HPP

#include "spaecies.hpp"

#include "rainshaft_integrator.hpp"

class FixedSubstepIntegrator : public RainshaftIntegrator {

public:

  // SPS: The usual issue here where we should use a different pointer type.
  FixedSubstepIntegrator(const RainshaftIntegrator *inner_integrator,
                         double dt_in);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const spaecies::State<double>& initial_state) const;

  const double dt;

private:

  const RainshaftIntegrator * integrator;

};

#endif // FIXED_SUBSTEP_INTEGRATOR_HPP
