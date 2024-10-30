#ifndef OPERATOR_SPLITTING_INTEGRATOR_HPP
#define OPERATOR_SPLITTING_INTEGRATOR_HPP

#include "sundials_integrator.hpp"

class OperatorSplittingIntegrator : public SundialsIntegrator<2>
{
private:
  const double dt;
  const double dt_exp;
  const double dt_imp;
  const int order;

public:
  OperatorSplittingIntegrator(const RainshaftConstants &constants,
                                       const RainshaftGrid &grid,
                                       const RainshaftProcess *const process_exp,
                                       const RainshaftProcess *const process_imp,
                                       const VarDescList& state_descs,
                                       const VarDescList& tend_descs,
                                       const double dt,
                                       const double dt_exp = 0,
                                       const double dt_imp = 0,
                                       const int order = 1,
                                       const int steps_per_output = -1);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const StateConst &initial_state) const;
};

#endif // OPERATOR_SPLITTING_INTEGRATOR_HPP
