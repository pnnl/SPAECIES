#ifndef OPERATOR_SPLITTING_INTEGRATOR_HPP
#define OPERATOR_SPLITTING_INTEGRATOR_HPP

#include "sundials_integrator.hpp"
#include "sedimentation.hpp"

class OperatorSplittingIntegrator : public SundialsIntegrator<2>
{
private:
  const double dt;
  const double dt_partition_1;
  const double dt_partition_2;
  const bool cfl_substep;
  const int order;
  const double rel_tol;
  const bool postprocess;

public:
  OperatorSplittingIntegrator(const RainshaftConstants &constants,
                              const RainshaftGrid &grid,
                              const SizeLimiters &size_limiters,
                              const Sedimentation *const process_partition_1,
                              const RainshaftProcess *const process_partition_2,
                              const VarDescList &state_descs,
                              const VarDescList &tend_descs,
                              const State& abs_tol,
                              const double dt,
                              const double dt_partition_1 = 0,
                              const double dt_partition_2 = 0,
                              const bool cfl_substep = false,
                              const int order = 1,
                              const double rel_tol = 1.e-4,
                              const bool postprocess = false,
                              const bool regularize_lambdar = true,
                              const int steps_per_output = -1);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const StateConst &initial_state,
                              int& error_flag) const;
};

#endif // OPERATOR_SPLITTING_INTEGRATOR_HPP
