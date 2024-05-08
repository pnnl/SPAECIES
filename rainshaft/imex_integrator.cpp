#include "imex_integrator.hpp"
#include "arkode/arkode_arkstep.h"
#include "nvector/nvector_serial.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

IMEXIntegrator::IMEXIntegrator(const RainshaftConstants &constants,
                                       const RainshaftGrid &grid,
                                       const RainshaftProcess *const process_exp,
                                       const RainshaftProcess *const process_imp,
                                       const double dt,
                                       const int order,
                                       const int steps_per_output)
    : SundialsIntegrator(constants, grid, {process_exp, process_imp}, steps_per_output),
      dt(dt), order(order)
{
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution IMEXIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const RainshaftState &initial_state) const
{
  auto y = state_to_n_vector(sun_ctxt, initial_state);
  void *arkode_mem = ARKStepCreate(create_f<0>(), create_f<1>(), initial_time, y, sun_ctxt);
  ARKStepSetUserData(arkode_mem, (void *)&user_data);
  ARKStepSetFixedStep(arkode_mem, dt);
  ARKStepSetOrder(arkode_mem, order);
  ARKStepSetMaxNumSteps(arkode_mem, -1); // Set no limit on the number of steps
  ARKStepSetStopTime(arkode_mem, final_time);

  SUNLinearSolver LS = nullptr; 
  SUNMatrix J = nullptr;
  J = SUNDenseMatrix(N_VGetLength(y), N_VGetLength(y), sun_ctxt);
  LS = SUNLinSol_Dense(y, J, sun_ctxt);
  ARKStepSetLinearSolver(arkode_mem, LS, J);
  // ARKStepSetMaxNonlinIters(arkode_mem, 160);
  // ARKStepSetNonlinConvCoef(arkode_mem, SUN_RCONST(1.e-3));
  ARKStepSetSafetyFactor(arkode_mem, 0.8);

  // ARKODE currently approximates the optimal step size with (err)^(-1/(embedded order)), but it's more common to use
  // (err)^(-1/(min(embedded order, order) + 1)). This adjusts the exponent to use the latter.
  ARKStepSetAdaptivityAdjustment(arkode_mem, 0);

  const sunrealtype fac = 1.;
  const sunrealtype reltol = fac * 1.e-2;
  auto abstol = N_VClone(y);
  auto tol_data = N_VGetArrayPointer_Serial(abstol);
  const auto nz = initial_state.t.size();
  for (sunindextype j = 0; j != nz; ++j)
  {
    tol_data[j] = fac * 1.e-1;
  }
  for (sunindextype j = 0; j != nz; ++j)
  {
    tol_data[nz + j] = fac * 1.e-5;
  }
  for (sunindextype j = 0; j != nz; ++j)
  {
    tol_data[2 * nz + j] = fac * 1.e-1;
  }
  for (sunindextype j = 0; j != nz; ++j)
  {
    tol_data[3 * nz + j] = fac * 1.e-8;
  }
  ARKStepSVtolerances(arkode_mem, reltol, abstol);
  N_VDestroy(abstol);

  auto solution = evolve(
      ARKStepEvolve,
      [arkode_mem]()
      {
        long int nfe = 0;
        long int nfi = 0;
        ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
        return nfe + nfi;
      },
      arkode_mem,
      final_time,
      y,
      ARK_NORMAL,
      ARK_ONE_STEP);

  N_VDestroy(y);
  // SPS: Make RAII wrapper for this.
  ARKStepFree(&arkode_mem);
  SUNMatDestroy(J);
  SUNLinSolFree(LS);
  return solution;
}
