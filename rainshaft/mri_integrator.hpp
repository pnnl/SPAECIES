#ifndef MRI_INTEGRATOR_HPP
#define MRI_INTEGRATOR_HPP
#include "sundials_integrator.hpp"

class MRIIntegrator : public SundialsIntegrator<3>
{

private:
  const double dt_fast;
  const double dt_slow;
  const int order;

public:
  MRIIntegrator(const RainshaftConstants &constants,
                const RainshaftGrid &grid,
                const RainshaftProcess *const process_fast,
                const RainshaftProcess *const process_slow_exp,
                const RainshaftProcess *const process_slow_imp,
                const std::vector<spaecies::VarDescPtr>& var_descs,
                const double dt_fast,
                const double dt_slow,
                const int order = 3,
                const int steps_per_output = -1);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const spaecies::VariableArray<double> &initial_state) const;
};

#endif // MRI_INTEGRATOR_HPP
