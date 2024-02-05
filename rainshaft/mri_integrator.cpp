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
                             const std::vector<spaecies::VarDescPtr>& var_descs,
                             const double dt_fast,
                             const double dt_slow,
                             const int order,
                             const int steps_per_output)
    : SundialsIntegrator(constants, grid, {process_fast, process_slow_exp, process_slow_imp}, var_descs, steps_per_output),
      dt_fast(dt_fast), dt_slow(dt_slow), order(order)
{
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution MRIIntegrator::integrate(double initial_time,
                                           double final_time,
                                           const spaecies::VariableArray<double> &initial_state) const
{
  auto y_init = state_to_n_vector(sun_ctxt, initial_state);

  /* create an ARKStep object, setting fast (inner) right-hand side
     functions and the initial condition */
  auto inner_arkode_mem = ARKStepCreate(create_f<0>(), nullptr, initial_time, y_init, sun_ctxt);
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
  void *outer_arkode_mem = MRIStepCreate(create_f<1>(), create_f<2>(), initial_time, y_init, stepper, sun_ctxt);

  /* Pass udata to user functions */
  ARKodeSetUserData(outer_arkode_mem, (void *)&user_data);
  ARKodeSetOrder(outer_arkode_mem, order);
  ARKodeSetMaxNumSteps(outer_arkode_mem, -1);
  ARKodeSetStopTime(outer_arkode_mem, final_time);

  /* Set the slow step size */
  ARKodeSetFixedStep(outer_arkode_mem, dt_slow);

  auto solution = evolve(
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
      outer_arkode_mem,
      final_time,
      y_init,
      ARK_NORMAL,
      ARK_ONE_STEP);

  N_VDestroy(y_init);
  MRIStepInnerStepper_Free(&stepper);
  ARKodeFree(&inner_arkode_mem);
  ARKodeFree(&outer_arkode_mem);

  return solution;
}
