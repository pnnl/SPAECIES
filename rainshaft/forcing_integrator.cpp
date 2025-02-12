#include "forcing_integrator.hpp"
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_forcingstep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_stepper.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_lapackdense.h"

ForcingIntegrator::ForcingIntegrator(const RainshaftConstants &constants,
                                     const RainshaftGrid &grid,
                                     const RainshaftProcess *const process_unforced,
                                     const RainshaftProcess *const process_forced,
                                     const VarDescList &state_descs,
                                     const VarDescList &tend_descs,
                                     const double dt,
                                     const double dt_unforced,
                                     const double dt_forced,
                                     const bool postprocess,
                                     const int steps_per_output)
    : SundialsIntegrator(constants, grid, {process_unforced, process_forced}, state_descs, tend_descs, steps_per_output),
      dt(dt), dt_unforced(dt_unforced), dt_forced(dt_forced), postprocess(postprocess)
{
  if (dt == 0.0 || dt_unforced == 0.0 || dt_forced == 0)
  {
    throw std::invalid_argument("The forcing method does not support adaptivity");
  }
}

RainshaftSolution ForcingIntegrator::integrate(double initial_time,
                                               double final_time,
                                               const StateConst &initial_state) const
{
  const N_Vector y = view_to_n_vector(sun_ctxt, initial_state);

  void *unforced_mem = ERKStepCreate(create_f<0>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(unforced_mem, (void *)&user_data);
  ARKodeSetFixedStep(unforced_mem, dt_unforced);
  ARKodeSetOrder(unforced_mem, 1);
  ARKodeSetMaxNumSteps(unforced_mem, -1); // Set no limit on the number of steps
  if (postprocess)
  {
    ARKodeSetPostprocessStageFn(unforced_mem, postprocess_positive);
    ARKodeSetPostprocessStepFn(unforced_mem, postprocess_positive);
  }

  SUNAdaptController unforced_controller = SUNAdaptController_I(sun_ctxt);
  ARKodeSetAdaptController(unforced_mem, unforced_controller);
  SUNAdaptController_SetErrorBias(unforced_controller, 1);
  ARKodeSetSafetyFactor(unforced_mem, 0.9);
  ARKodeSetAdaptivityAdjustment(unforced_mem, 0);
  ARKodeSetFixedStepBounds(unforced_mem, 1, 1);

  void *forced_mem = ERKStepCreate(create_f<1>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(forced_mem, (void *)&user_data);
  ARKodeSetFixedStep(forced_mem, dt_forced);
  ARKodeSetOrder(forced_mem, 1);
  ARKodeSetMaxNumSteps(forced_mem, -1); // Set no limit on the number of steps
  if (postprocess)
  {
    ARKodeSetPostprocessStageFn(forced_mem, postprocess_positive);
    ARKodeSetPostprocessStepFn(forced_mem, postprocess_positive);
  }

  SUNAdaptController forced_controller = SUNAdaptController_I(sun_ctxt);
  ARKodeSetAdaptController(forced_mem, forced_controller);
  SUNAdaptController_SetErrorBias(forced_controller, 1);
  ARKodeSetSafetyFactor(forced_mem, 0.9);
  ARKodeSetAdaptivityAdjustment(forced_mem, 0);
  ARKodeSetFixedStepBounds(forced_mem, 1, 1);

  std::array<SUNStepper, 2> steppers{};
  ARKodeCreateSUNStepper(unforced_mem, &steppers[0]);
  ARKodeCreateSUNStepper(forced_mem, &steppers[1]);
  void *splitting_mem = ForcingStepCreate(steppers[0], steppers[1], initial_time, y, sun_ctxt);
  ARKodeSetFixedStep(splitting_mem, dt);
  ARKodeSetMaxNumSteps(splitting_mem, -1); // Set no limit on the number of steps
  ARKodeSetStopTime(splitting_mem, final_time);

  const RainshaftSolution solution = evolve(
      ARKodeEvolve,
      [unforced_mem, forced_mem]()
      {
        long int num_unforced_rhs_evals = 0;
        ARKodeGetNumRhsEvals(unforced_mem, 0, &num_unforced_rhs_evals);
        long int num_forced_rhs_evals = 0;
        ARKodeGetNumRhsEvals(forced_mem, 0, &num_forced_rhs_evals);
        return num_unforced_rhs_evals + num_forced_rhs_evals;
      },
      []() {},
      splitting_mem,
      final_time,
      y,
      ARK_NORMAL,
      ARK_ONE_STEP);

  N_VDestroy(y);
  SUNAdaptController_Destroy(unforced_controller);
  SUNAdaptController_Destroy(forced_controller);
  SUNStepper_Destroy(&steppers[0]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&unforced_mem);
  ARKodeFree(&forced_mem);
  ARKodeFree(&splitting_mem);
  return solution;
}
