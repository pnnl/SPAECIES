#include "mri_integrator.hpp"
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_mristep.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include "sundials/sundials_context.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#include "sedimentation.hpp"
#include <algorithm>
#include <cmath>

MRIIntegrator::MRIIntegrator(const RainshaftConstants &constants,
                             const RainshaftGrid &grid,
                             const RainshaftProcess *const process_fast,
                             const RainshaftProcess *const process_slow_exp,
                             const RainshaftProcess *const process_slow_imp,
                             const VarDescList& state_descs,
                             const VarDescList& tend_descs,
                             const double dt_fast,
                             const double dt_slow,
                             const int order,
                             const double rel_tol,
                             const int steps_per_output)
    : SundialsIntegrator(constants, grid, {process_fast, process_slow_exp, process_slow_imp}, state_descs, tend_descs, steps_per_output),
      dt_fast(dt_fast), dt_slow(dt_slow), order(order), rel_tol(rel_tol)
{
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution MRIIntegrator::integrate(double initial_time,
                                           double final_time,
                                           const StateConst &initial_state) const
{
  const N_Vector y = view_to_n_vector(sun_ctxt, initial_state);

  /* create an ARKStep object, setting fast (inner) right-hand side
     functions and the initial condition */
  void *inner_arkode_mem = ARKStepCreate(create_f<0>(), nullptr, initial_time, y, sun_ctxt);
  ARKodeSetOrder(inner_arkode_mem, order);
  ARKodeSetUserData(inner_arkode_mem, (void *)&user_data);
  ARKodeSetMaxNumSteps(inner_arkode_mem, -1);

  // fixed step for fast solver
  ARKodeSetFixedStep(inner_arkode_mem, dt_fast);

  /* create MRIStepInnerStepper wrapper for the ARKStep memory block */
  MRIStepInnerStepper stepper = nullptr;
  ARKStepCreateMRIStepInnerStepper(inner_arkode_mem, &stepper);

  /* create an MRIStep object, setting the slow (outer) right-hand side
     functions and the initial condition */
  void *outer_arkode_mem = MRIStepCreate(create_f<1>(), create_f<2>(), initial_time, y, stepper, sun_ctxt);

  /* Pass udata to user functions */
  ARKodeSetUserData(outer_arkode_mem, (void *)&user_data);
  ARKodeSetOrder(outer_arkode_mem, order);
  ARKodeSetMaxNumSteps(outer_arkode_mem, -1);
  ARKodeSetStopTime(outer_arkode_mem, final_time);

  /* Set the slow step size */
  ARKodeSetFixedStep(outer_arkode_mem, dt_slow);

  const RainshaftSolution solution = evolve(
      ARKodeEvolve,
      [outer_arkode_mem, inner_arkode_mem]()
      {
        long int nfse = 0;
        long int nfsi = 0;
        MRIStepGetNumRhsEvals(outer_arkode_mem, &nfse, &nfsi);

        long int nffe = 0;
        long int nffi = 0;
        ARKStepGetNumRhsEvals(inner_arkode_mem, &nffe, &nffi);
        return nfse + nfsi + nffe + nffi;
      },
      []() {},
      outer_arkode_mem,
      final_time,
      y,
      ARK_NORMAL,
      ARK_ONE_STEP);

  N_VDestroy(y);
  MRIStepInnerStepper_Free(&stepper);
  ARKodeFree(&inner_arkode_mem);
  ARKodeFree(&outer_arkode_mem);

  return solution;
}
