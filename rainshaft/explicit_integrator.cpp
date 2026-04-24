#include "explicit_integrator.hpp"
#include "arkode/arkode_erkstep.h"
#include "nvector/nvector_serial.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"

ExplicitIntegrator::ExplicitIntegrator(const RainshaftConstants &constants,
                                       const RainshaftGrid &grid,
                                       const SizeLimiters &size_limiters,
                                       const RainshaftProcess *const process,
                                       const VarDescList& state_descs,
                                       const VarDescList& tend_descs,
                                       const State& abs_tol,
                                       const double dt,
                                       const int order,
                                       const double rel_tol,
                                       const bool postprocess,
                                       const bool regularize_lambdar,
                                       const int steps_per_output)
    : SundialsIntegrator(constants, grid, size_limiters, {process}, state_descs, tend_descs, abs_tol, steps_per_output, regularize_lambdar),
      dt(dt), order(order), rel_tol(rel_tol), postprocess(postprocess)
{
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution ExplicitIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const StateConst& initial_state,
                                                int& error_flag) const
{
  const N_Vector y = view_to_n_vector(sun_ctxt, initial_state);
  void *arkode_mem = ERKStepCreate(create_f<0>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(arkode_mem, (void *)&user_data);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetOrder(arkode_mem, order);
  ARKodeSetMaxNumSteps(arkode_mem, -1); // Set no limit on the number of steps
  ARKodeSetStopTime(arkode_mem, final_time);
  if (postprocess) {
    ARKodeSetPostprocessStageFn(arkode_mem, postprocess_positive);
    ARKodeSetPostprocessStepFn(arkode_mem, postprocess_positive);
  }

  const N_Vector abs_tol = create_abs_tol_n_vector();
  ARKodeSVtolerances(arkode_mem, rel_tol, abs_tol);
  N_VDestroy(abs_tol);

  const RainshaftSolution solution = evolve(
      ARKodeEvolve,
      [arkode_mem]()
      {
        long int num_rhs_evals = 0;
        ARKodeGetNumRhsEvals(arkode_mem, 0, &num_rhs_evals);
        return num_rhs_evals;
      },
      []() {},
      arkode_mem,
      final_time,
      y,
      ARK_NORMAL,
      ARK_ONE_STEP,
      error_flag);

  N_VDestroy(y);
  // SPS: Make RAII wrapper for this.
  ARKodeFree(&arkode_mem);
  return solution;
}
