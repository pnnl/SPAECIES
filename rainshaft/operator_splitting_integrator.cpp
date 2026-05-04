#include "operator_splitting_integrator.hpp"
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_splittingstep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_stepper.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_lapackdense.h"

OperatorSplittingIntegrator::OperatorSplittingIntegrator(const RainshaftConstants &constants,
                                                         const RainshaftGrid &grid,
                                                         const SizeLimiters &size_limiters,
                                                         const Sedimentation *const process_partition_1,
                                                         const RainshaftProcess *const process_partition_2,
                                                         const VarDescList &state_descs,
                                                         const VarDescList &tend_descs,
                                                         const State& abs_tol,
                                                         const double dt,
                                                         const double dt_partition_1,
                                                         const double dt_partition_2,
                                                         const bool cfl_substep,
                                                         const int order,
                                                         const double rel_tol,
                                                         const bool postprocess,
                                                         const bool regularize_lambdar,
                                                         const int steps_per_output)
    : SundialsIntegrator(constants, grid, size_limiters, {process_partition_1, process_partition_2}, state_descs, tend_descs, abs_tol, steps_per_output, regularize_lambdar),
      dt(dt), dt_partition_1(dt_partition_1), dt_partition_2(dt_partition_2), cfl_substep(cfl_substep), order(order), rel_tol(rel_tol), postprocess(postprocess)
{
}

static void setPartitionOrder(void *partition_mem, int order, bool cfl_substep=false) {
  if (order == 1 && cfl_substep) {
      /* A bit of a hack to get 1st order methods working. Currently ARKODE
       * requires the method to have an embedding to use a stability function */
      sunrealtype a[] = {0.0};
      sunrealtype b[] = {1.0};
      ARKodeButcherTable euler = ARKodeButcherTable_Create(1, 1, 1, a, a, b, b);
      ERKStepSetTable(partition_mem, euler);
      ARKodeButcherTable_Free(euler);
  } else if (order == 2) {
    /* The default 2nd and 3rd order ERK methods have the first same as last
     * (FSAL) property. This can't be exploited with operator splitting so they
     * effectively have a redundant stage. We override the defaults with methods
     * without the FSAL property. */
    ERKStepSetTableNum(partition_mem, ARKODE_RALSTON_EULER_2_1_2);
  } else if (order == 3) {
    ERKStepSetTableNum(partition_mem, ARKODE_SHU_OSHER_3_2_3);
  } else {
    ARKodeSetOrder(partition_mem, order);
  }
}

RainshaftSolution OperatorSplittingIntegrator::integrate(double initial_time,
                                                         double final_time,
                                                         const StateConst &initial_state,
                                                         int& error_flag) const
{
  const N_Vector y = view_to_n_vector(sun_ctxt, initial_state);

  void *partition_1_mem = ERKStepCreate(create_f<0>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(partition_1_mem, (void *)&user_data);
  ARKodeSetFixedStep(partition_1_mem, dt_partition_1);
  setPartitionOrder(partition_1_mem, order, cfl_substep);
  ARKodeSetMaxNumSteps(partition_1_mem, -1); // Set no limit on the number of steps
  if (cfl_substep)
  {
    /* Currently a CFL fraction of 1 isn't allowed so use the largest floating
     * point number less than 1. This should be resolved in the next SUNDIALS
     * release */
    ARKodeSetCFLFraction(partition_1_mem, std::nextafter(1.0, 0.0));
    ARKodeSetStabilityFn(partition_1_mem, create_stability<0, Sedimentation>(), (void*)&user_data);
  }
  if (postprocess)
  {
    ARKodeSetPostprocessStageFn(partition_1_mem, postprocess_positive);
    ARKodeSetPostprocessStepFn(partition_1_mem, postprocess_positive);
  }

  void *partition_2_mem = ERKStepCreate(create_f<1>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(partition_2_mem, (void *)&user_data);
  ARKodeSetFixedStep(partition_2_mem, dt_partition_2);
  setPartitionOrder(partition_2_mem, order);
  ARKodeSetMaxNumSteps(partition_2_mem, -1); // Set no limit on the number of steps
  if (postprocess)
  {
    ARKodeSetPostprocessStageFn(partition_2_mem, postprocess_positive);
    ARKodeSetPostprocessStepFn(partition_2_mem, postprocess_positive);
  }

  const N_Vector abs_tol = create_abs_tol_n_vector();
  ARKodeSVtolerances(partition_1_mem, rel_tol, abs_tol);
  ARKodeSVtolerances(partition_2_mem, rel_tol, abs_tol);
  N_VDestroy(abs_tol);

  std::array<SUNStepper, 2> steppers{};
  ARKodeCreateSUNStepper(partition_1_mem, &steppers[1]);
  ARKodeCreateSUNStepper(partition_2_mem, &steppers[0]);
  void *splitting_mem = SplittingStepCreate(steppers.data(), steppers.size(), initial_time, y, sun_ctxt);
  ARKodeSetFixedStep(splitting_mem, dt);
  ARKodeSetOrder(splitting_mem, order);
  ARKodeSetMaxNumSteps(splitting_mem, -1); // Set no limit on the number of steps
  ARKodeSetStopTime(splitting_mem, final_time);

  const RainshaftSolution solution = evolve(
      ARKodeEvolve,
      [partition_1_mem, partition_2_mem]()
      {
        long int num_partition_1_evals = 0;
        ARKodeGetNumRhsEvals(partition_1_mem, 0, &num_partition_1_evals);
        long int num_partition_2_evals = 0;
        ARKodeGetNumRhsEvals(partition_2_mem, 0, &num_partition_2_evals);
        return num_partition_1_evals + num_partition_2_evals;
      },
      []() {},
      splitting_mem,
      final_time,
      y,
      ARK_NORMAL,
      ARK_ONE_STEP,
      error_flag);

  N_VDestroy(y);
  SUNStepper_Destroy(&steppers[0]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&partition_1_mem);
  ARKodeFree(&partition_2_mem);
  ARKodeFree(&splitting_mem);
  return solution;
}
