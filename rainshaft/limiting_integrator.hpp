#ifndef LIMITING_INTEGRATOR_HPP
#define LIMITING_INTEGRATOR_HPP

#include "rainshaft_integrator.hpp"
#include "size_limiters.hpp"

class LimitingIntegrator : public RainshaftIntegrator {

public:

  LimitingIntegrator(const RainshaftConstants &constants,
                     const SizeLimiters& size_limiters,
                     const RainshaftIntegrator& sub_integrator);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const StateConst& initial_state,
                              int& error_flag) const;

private:

  const RainshaftConstants &constants;

  const SizeLimiters size_limiters;

  const RainshaftIntegrator &integrator;

};

#endif // LIMITING_INTEGRATOR_HPP
