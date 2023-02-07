#ifndef RAINSHAFT_INTEGRATOR_HPP
#define RAINSHAFT_INTEGRATOR_HPP
#include "rainshaft_state.hpp"
#include "rainshaft_tendency.hpp"

class RainshaftIntegrator {

public:

  // Constructor requires timestep.
  RainshaftIntegrator(double dt_in);

  RainshaftState apply_tend(const RainshaftState& state,
                            const RainshaftTendency& tend);

  const double dt;

};

#endif // RAINSHAFT_INTEGRATOR_HPP
