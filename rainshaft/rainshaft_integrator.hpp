#ifndef RAINSHAFT_INTEGRATOR_HPP
#define RAINSHAFT_INTEGRATOR_HPP

#include "rainshaft_solution.hpp"
#include "rainshaft_types.hpp"

class RainshaftIntegrator {

public:

  virtual RainshaftSolution integrate(double initial_time,
                                      double final_time,
                                      const StateConst& initial_state) const = 0;

  virtual ~RainshaftIntegrator() {}

};

#endif // RAINSHAFT_INTEGRATOR_HPP
