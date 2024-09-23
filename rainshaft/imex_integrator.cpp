#include "imex_integrator.hpp"
#include "arkode/arkode_arkstep.h"
#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_lapackdense.h"

IMEXIntegrator::IMEXIntegrator(const RainshaftConstants &constants,
                               const RainshaftGrid &grid,
                               const RainshaftProcess *const process_exp,
                               const RainshaftProcess *const process_imp,
                                       const VarDescList& state_descs,
                                       const VarDescList& tend_descs,
                               const double dt,
                               const int order,
                               const int steps_per_output,
                               const std::optional<std::string> jacobian_file)
    : SundialsIntegrator(constants, grid, {process_exp, process_imp}, state_descs, tend_descs, steps_per_output),
      dt(dt), order(order), jacobian_file(jacobian_file)
{
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution IMEXIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const StateConst &initial_state) const
{
  N_Vector y = state_to_n_vector(sun_ctxt, initial_state);
  void *arkode_mem = ARKStepCreate(create_f<0>(), create_f<1>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(arkode_mem, (void *)&user_data);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetOrder(arkode_mem, order);
  ARKodeSetMaxNumSteps(arkode_mem, -1); // Set no limit on the number of steps
  ARKodeSetStopTime(arkode_mem, final_time);

  SUNMatrix jac = SUNDenseMatrix(N_VGetLength(y), N_VGetLength(y), sun_ctxt);
  SUNLinearSolver LS = SUNLinSol_LapackDense(y, jac, sun_ctxt);
  ARKodeSetLinearSolver(arkode_mem, LS, jac);
  ARKodeSetJacFn(arkode_mem, create_jac<1>());

  auto controller = SUNAdaptController_I(sun_ctxt);
  ARKodeSetAdaptController(arkode_mem, controller);
  // Use no bias and instead rely on a safety factor. Shampine uses these values in the MATLAB ODE suite
  SUNAdaptController_SetErrorBias(controller, 1);
  ARKodeSetSafetyFactor(arkode_mem, 0.9);
  // ARKODE currently approximates the optimal step size with (err)^(-1/(embedded order)), but it's more common to use
  // (err)^(-1/(min(embedded order, order) + 1)). This adjusts the exponent to use the latter.
  ARKodeSetAdaptivityAdjustment(arkode_mem, 0);
  ARKodeSetFixedStepBounds(arkode_mem, 1, 1); // Remove deadzone

  const sunrealtype fac = 1.;
  const sunrealtype reltol = fac * 1.e-2;
  auto abstol = N_VClone(y);
  auto tol_data = N_VGetArrayPointer_Serial(abstol);
  const auto nz = user_data.grid.nlev;
  for (std::size_t j = 0; j != nz; ++j)
  {
    tol_data[j] = fac * 1.e-1;
  }
  for (std::size_t j = 0; j != nz; ++j)
  {
    tol_data[nz + j] = fac * 1.e-5;
  }
  for (std::size_t j = 0; j != nz; ++j)
  {
    tol_data[2 * nz + j] = fac * 1.e-1;
  }
  for (std::size_t j = 0; j != nz; ++j)
  {
    tol_data[3 * nz + j] = fac * 1.e-8;
  }
  ARKodeSVtolerances(arkode_mem, reltol, abstol);
  N_VDestroy(abstol);

  auto solution = evolve(
      ARKodeEvolve,
      [arkode_mem]()
      {
        long int nfe = 0;
        long int nfi = 0;
        ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
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
        const auto filename = jacobian_file.value() + "." + std::to_string(steps);
        const auto file = fopen(filename.c_str(), "w");
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
