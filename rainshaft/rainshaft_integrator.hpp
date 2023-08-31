#ifndef RAINSHAFT_INTEGRATOR_HPP
#define RAINSHAFT_INTEGRATOR_HPP
#include "rainshaft_state.hpp"
#include "rainshaft_solution.hpp"

class RainshaftIntegrator {

public:

  virtual RainshaftSolution integrate(double initial_time,
                                      double final_time,
                                      const RainshaftState& initial_state) = 0;

};

#endif // RAINSHAFT_INTEGRATOR_HPP
