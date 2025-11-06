#include "imex_integrator.hpp"
#include "arkode/arkode_arkstep.h"
#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_lapackdense.h"

IMEXIntegrator::IMEXIntegrator(const RainshaftConstants &constants,
                               const RainshaftGrid &grid,
                               const SizeLimiters &size_limiters,
                               const RainshaftProcess *const process_imp,
                               const RainshaftProcess *const process_exp,
                               const VarDescList& state_descs,
                               const VarDescList& tend_descs,
                               const double dt,
                               const int order,
                               const double rel_tol,
                               const bool postprocess,
                               const bool regularize_lambdar,
                               const int steps_per_output,
                               const std::optional<std::string> jacobian_file)
    : SundialsIntegrator(constants, grid, size_limiters, {process_exp, process_imp}, state_descs, tend_descs, steps_per_output, regularize_lambdar),
      dt(dt), order(order), rel_tol(rel_tol), postprocess(postprocess), jacobian_file(jacobian_file)
{
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution IMEXIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const StateConst &initial_state,
                                                int& error_flag) const
{
  const N_Vector y = view_to_n_vector(sun_ctxt, initial_state);
  void *arkode_mem = ARKStepCreate(create_f<0>(), create_f<1>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(arkode_mem, (void *)&user_data);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetOrder(arkode_mem, order);
  ARKodeSetMaxNumSteps(arkode_mem, -1); // Set no limit on the number of steps
  ARKodeSetStopTime(arkode_mem, final_time);
  if (postprocess) {
    ARKodeSetPostprocessStageFn(arkode_mem, postprocess_positive);
    ARKodeSetPostprocessStepFn(arkode_mem, postprocess_positive);
  }

  SUNMatrix jac = SUNDenseMatrix(N_VGetLength(y), N_VGetLength(y), sun_ctxt);
  SUNLinearSolver ls = SUNLinSol_LapackDense(y, jac, sun_ctxt);
  ARKodeSetLinearSolver(arkode_mem, ls, jac);
  ARKodeSetJacFn(arkode_mem, create_jac<1>());
  ARKodeSetDeduceImplicitRhs(arkode_mem, true);

  if (dt > 0) {
    ARKodeSetMaxNonlinIters(arkode_mem, 1000);
  } else {
    ARKodeSetPredictorMethod(arkode_mem, 1);
  }

  const N_Vector abs_tol = fill_abs_tol_vector(N_VClone(y));
  ARKodeSVtolerances(arkode_mem, rel_tol, abs_tol);
  N_VDestroy(abs_tol);

  const RainshaftSolution solution = evolve(
      ARKodeEvolve,
      [arkode_mem]()
      {
        long int nfe = 0;
        long int nfi = 0;
        ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
        ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
        return nfe + nfi;
      },
      [&]()
      {
        if (!jacobian_file.has_value())
        {
          return;
        }
        SUNMatrix mat = nullptr;
        ARKodeGetJac(arkode_mem, &mat);
        if (mat == nullptr) {
          return;
        }
        
        long int steps = 0;
        ARKodeGetNumSteps(arkode_mem, &steps);
        const std::string filename = jacobian_file.value() + "." + std::to_string(steps);
        std::FILE* const file = fopen(filename.c_str(), "w");
        SUNDenseMatrix_Print(mat, file);
        fclose(file);
      },
      arkode_mem,
      final_time,
      y,
      ARK_NORMAL,
      ARK_ONE_STEP,
      error_flag);

  N_VDestroy(y);
  // SPS: Make RAII wrapper for this.
  ARKodeFree(&arkode_mem);
  SUNMatDestroy(jac);
  SUNLinSolFree(ls);
  return solution;
}
