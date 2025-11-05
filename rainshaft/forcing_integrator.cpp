#include "forcing_integrator.hpp"
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_forcingstep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_stepper.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_lapackdense.h"

ForcingIntegrator::ForcingIntegrator(const RainshaftConstants &constants,
                                     const RainshaftGrid &grid,
                                     const Sedimentation *const process_forced,
                                     const RainshaftProcess *const process_unforced,
                                     const VarDescList &state_descs,
                                     const VarDescList &tend_descs,
                                     const double dt,
                                     const double dt_forced,
                                     const double dt_unforced,
                                     const bool cfl_substep,
                                     const bool postprocess,
                                     const bool regularize_lambdar,
                                     const int steps_per_output)
    : SundialsIntegrator(constants, grid, {process_forced, process_unforced}, state_descs, tend_descs, steps_per_output, regularize_lambdar),
      dt(dt), dt_forced(dt_forced), dt_unforced(dt_unforced), cfl_substep(cfl_substep), postprocess(postprocess)
{
  if (dt == 0.0 || (dt_forced == 0.0 && !cfl_substep) || dt_unforced == 0)
  {
    throw std::invalid_argument("The forcing method does not support adaptivity");
  }
}

RainshaftSolution ForcingIntegrator::integrate(double initial_time,
                                               double final_time,
                                               const StateConst &initial_state,
                                               int& error_flag) const
{
  const N_Vector y = view_to_n_vector(sun_ctxt, initial_state);

  void *forced_mem = ERKStepCreate(create_f<0>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(forced_mem, (void *)&user_data);
  ARKodeSetFixedStep(forced_mem, dt_forced);

  if (cfl_substep) {
      /* Currently a CFL fraction of 1 isn't allowed so use the largest floating
      * point number less than 1. This should be resolved in the next SUNDIALS
      * release */
      ARKodeSetCFLFraction(forced_mem, std::nextafter(1.0, 0.0));
      ARKodeSetStabilityFn(forced_mem, create_stability<0, Sedimentation>(), (void*)&user_data);

      /* A bit of a hack to get 1st order methods working. Currently ARKODE
       * requires the method to have an embedding to use a stability function */
      sunrealtype a[] = {0.0};
      sunrealtype b[] = {1.0};
      ARKodeButcherTable euler = ARKodeButcherTable_Create(1, 1, 1, a, a, b, b);
      ERKStepSetTable(forced_mem, euler);
      ARKodeButcherTable_Free(euler);
  } else {
    ARKodeSetOrder(forced_mem, 1);
  }
  
  ARKodeSetMaxNumSteps(forced_mem, -1); // Set no limit on the number of steps
  if (postprocess)
  {
    ARKodeSetPostprocessStageFn(forced_mem, postprocess_positive);
    ARKodeSetPostprocessStepFn(forced_mem, postprocess_positive);
  }

  void *unforced_mem = ERKStepCreate(create_f<1>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(unforced_mem, (void *)&user_data);
  ARKodeSetFixedStep(unforced_mem, dt_unforced);
  ARKodeSetOrder(unforced_mem, 1);
  ARKodeSetMaxNumSteps(unforced_mem, -1); // Set no limit on the number of steps
  if (postprocess)
  {
    ARKodeSetPostprocessStageFn(unforced_mem, postprocess_positive);
    ARKodeSetPostprocessStepFn(unforced_mem, postprocess_positive);
  }

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
      ARK_ONE_STEP,
      error_flag);

  N_VDestroy(y);
  SUNStepper_Destroy(&steppers[0]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&unforced_mem);
  ARKodeFree(&forced_mem);
  ARKodeFree(&splitting_mem);
  return solution;
}
