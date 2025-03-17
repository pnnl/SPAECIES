#include "imex_integrator.hpp"
#include "arkode/arkode_arkstep.h"
#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_lapackdense.h"

IMEXIntegrator::IMEXIntegrator(const RainshaftConstants &constants,
                               const RainshaftGrid &grid,
                               const RainshaftProcess *const process_imp,
                               const RainshaftProcess *const process_exp,
                                       const VarDescList& state_descs,
                                       const VarDescList& tend_descs,
                               const double dt,
                               const int order,
                               const double rel_tol,
                               const int steps_per_output,
                               const bool postprocess,
                               const std::optional<std::string> jacobian_file)
    : SundialsIntegrator(constants, grid, {process_exp, process_imp}, state_descs, tend_descs, steps_per_output),
      dt(dt), order(order), rel_tol(rel_tol), postprocess(postprocess), jacobian_file(jacobian_file)
{
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution IMEXIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const StateConst &initial_state) const
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
  SUNLinearSolver LS = SUNLinSol_LapackDense(y, jac, sun_ctxt);
  ARKodeSetLinearSolver(arkode_mem, LS, jac);
  ARKodeSetJacFn(arkode_mem, create_jac<1>());
  ARKodeSetMaxNonlinIters(arkode_mem, 100);
  ARKodeSetDeduceImplicitRhs(arkode_mem, true);

  SUNAdaptController controller = SUNAdaptController_I(sun_ctxt);
  ARKodeSetAdaptController(arkode_mem, controller);
  // Use no bias and instead rely on a safety factor. Shampine uses these values in the MATLAB ODE suite
  SUNAdaptController_SetErrorBias(controller, 1);
  ARKodeSetSafetyFactor(arkode_mem, 0.9);
  // ARKODE currently approximates the optimal step size with (err)^(-1/(embedded order)), but it's more common to use
  // (err)^(-1/(min(embedded order, order) + 1)). This adjusts the exponent to use the latter.
  ARKodeSetAdaptivityAdjustment(arkode_mem, 0);
  ARKodeSetFixedStepBounds(arkode_mem, 1, 1); // Remove deadzone

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
      ARK_ONE_STEP);

  N_VDestroy(y);
  SUNAdaptController_Destroy(controller);
  // SPS: Make RAII wrapper for this.
  ARKodeFree(&arkode_mem);
  SUNMatDestroy(jac);
  SUNLinSolFree(LS);
  return solution;
}
