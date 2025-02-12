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
                                                         const RainshaftProcess *const process_partition_1,
                                                         const RainshaftProcess *const process_partition_2,
                                                         const VarDescList &state_descs,
                                                         const VarDescList &tend_descs,
                                                         const double dt,
                                                         const double dt_partition_1,
                                                         const double dt_partition_2,
                                                         const int order,
                                                         const double rel_tol,
                                                         const bool postprocess,
                                                         const int steps_per_output)
    : SundialsIntegrator(constants, grid, {process_partition_1, process_partition_2}, state_descs, tend_descs, steps_per_output),
      dt(dt), dt_partition_1(dt_partition_1), dt_partition_2(dt_partition_2), order(order), rel_tol(rel_tol), postprocess(postprocess)
{
}

RainshaftSolution OperatorSplittingIntegrator::integrate(double initial_time,
                                                         double final_time,
                                                         const StateConst &initial_state) const
{
  const N_Vector y = view_to_n_vector(sun_ctxt, initial_state);

  void *partition_1_mem = ERKStepCreate(create_f<0>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(partition_1_mem, (void *)&user_data);
  ARKodeSetFixedStep(partition_1_mem, dt_partition_1);
  ARKodeSetOrder(partition_1_mem, order);
  ARKodeSetMaxNumSteps(partition_1_mem, -1); // Set no limit on the number of steps
  if (postprocess)
  {
    ARKodeSetPostprocessStageFn(partition_1_mem, postprocess_positive);
    ARKodeSetPostprocessStepFn(partition_1_mem, postprocess_positive);
  }

  SUNAdaptController partition_1_controller = SUNAdaptController_I(sun_ctxt);
  ARKodeSetAdaptController(partition_1_mem, partition_1_controller);
  SUNAdaptController_SetErrorBias(partition_1_controller, 1);
  ARKodeSetSafetyFactor(partition_1_mem, 0.9);
  ARKodeSetAdaptivityAdjustment(partition_1_mem, 0);
  ARKodeSetFixedStepBounds(partition_1_mem, 1, 1);

  void *partition_2_mem = ERKStepCreate(create_f<1>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(partition_2_mem, (void *)&user_data);
  ARKodeSetFixedStep(partition_2_mem, dt_partition_2);
  ARKodeSetOrder(partition_2_mem, order);
  ARKodeSetMaxNumSteps(partition_2_mem, -1); // Set no limit on the number of steps
  if (postprocess)
  {
    ARKodeSetPostprocessStageFn(partition_2_mem, postprocess_positive);
    ARKodeSetPostprocessStepFn(partition_2_mem, postprocess_positive);
  }

  SUNAdaptController partition_2_controller = SUNAdaptController_I(sun_ctxt);
  ARKodeSetAdaptController(partition_2_mem, partition_2_controller);
  SUNAdaptController_SetErrorBias(partition_2_controller, 1);
  ARKodeSetSafetyFactor(partition_2_mem, 0.9);
  ARKodeSetAdaptivityAdjustment(partition_2_mem, 0);
  ARKodeSetFixedStepBounds(partition_2_mem, 1, 1);

  const N_Vector abs_tol = fill_abs_tol_vector(N_VClone(y));
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
      ARK_ONE_STEP);

  N_VDestroy(y);
  SUNAdaptController_Destroy(partition_1_controller);
  SUNAdaptController_Destroy(partition_2_controller);
  SUNStepper_Destroy(&steppers[0]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&partition_1_mem);
  ARKodeFree(&partition_2_mem);
  ARKodeFree(&splitting_mem);
  return solution;
}
