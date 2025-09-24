#ifndef SEQUENTIAL_SPLIT_INTEGRATOR_HPP
#define SEQUENTIAL_SPLIT_INTEGRATOR_HPP

#include "rainshaft_integrator.hpp"

class SequentialSplitIntegrator : public RainshaftIntegrator {

public:

  SequentialSplitIntegrator(const std::vector<const RainshaftIntegrator *>& sub_integrators);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const StateConst& initial_state,
                              int& error_flag) const;

private:

  const std::vector<const RainshaftIntegrator *> integrators;

};

#endif // SEQUENTIAL_SPLIT_INTEGRATOR_HPP
