#ifndef RAINSHAFT_INTEGRATOR_HPP
#define RAINSHAFT_INTEGRATOR_HPP

#include "spaecies.hpp"

#include "rainshaft_solution.hpp"

class RainshaftIntegrator {

public:

  virtual RainshaftSolution integrate(double initial_time,
                                      double final_time,
                                      const spaecies::State<const double>& initial_state) const = 0;

};

#endif // RAINSHAFT_INTEGRATOR_HPP
