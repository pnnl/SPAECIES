#ifndef EXPLICIT_INTEGRATOR_HPP
#define EXPLICIT_INTEGRATOR_HPP

#include "spaecies.hpp"

#include "sundials_integrator.hpp"

class ExplicitIntegrator : public SundialsIntegrator<1>
{
private:
  const double dt;
  const int order;

public:
  ExplicitIntegrator(const RainshaftConstants &constants,
                     const RainshaftGrid &grid,
                     const RainshaftProcess *const process,
                     const std::vector<spaecies::VarDescPtr>& var_descs,
                     const double dt = 0,
                     const int order = 4,
                     const int steps_per_output = -1);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const spaecies::VariableArray<double> &initial_state) const;
};

#endif // EXPLICIT_INTEGRATOR_HPP
