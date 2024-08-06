#include "imex_integrator.hpp"
#include "arkode/arkode_arkstep.h"
#include "nvector/nvector_serial.h"
#include "sunlinsol/sunlinsol_spbcgs.h"
#include "sunlinsol/sunlinsol_sptfqmr.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_lapackdense.h"
#include "sunlinsol/sunlinsol_band.h"
#include "sunmatrix/sunmatrix_band.h"
#include "sunlinsol/sunlinsol_lapackband.h"
#include "sunlinsol/sunlinsol_spgmr.h"

#include <iostream>

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
  ARKodeSetUserData(arkode_mem, (void *)&user_data);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetOrder(arkode_mem, order);
  ARKodeSetMaxNumSteps(arkode_mem, -1); // Set no limit on the number of steps
  ARKodeSetStopTime(arkode_mem, final_time);

  // SUNLinearSolver LS = SUNLinSol_SPTFQMR(y, SUN_PREC_NONE, 40, sun_ctxt);
  // ARKodeSetLinearSolver(arkode_mem, LS, nullptr);
  // ARKodeSetJacTimes(arkode_mem, nullptr, create_jac_prod<1>());

  // SUNLinearSolver LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 40, sun_ctxt);
  // ARKodeSetLinearSolver(arkode_mem, LS, nullptr);
  // ARKodeSetJacTimes(arkode_mem, nullptr, create_jac_prod<1>());

  // SUNMatrix jac = SUNDenseMatrix(N_VGetLength(y), N_VGetLength(y), sun_ctxt);
  // SUNLinearSolver LS = SUNLinSol_Dense(y, jac, sun_ctxt);
  // ARKodeSetLinearSolver(arkode_mem, LS, jac);
  // ARKodeSetJacFn(arkode_mem, create_jac<1>());

  SUNMatrix jac = SUNDenseMatrix(N_VGetLength(y), N_VGetLength(y), sun_ctxt);
  SUNLinearSolver LS = SUNLinSol_LapackDense(y, jac, sun_ctxt);
  ARKodeSetLinearSolver(arkode_mem, LS, jac);
  ARKodeSetJacFn(arkode_mem, create_jac<1>());

  // SUNMatrix jac = SUNBandMatrix(N_VGetLength(y), 0, 1, sun_ctxt);
  // SUNLinearSolver LS = SUNLinSol_Band(y, jac, sun_ctxt);
  // ARKodeSetLinearSolver(arkode_mem, LS, jac);
  // ARKodeSetJacFn(arkode_mem, create_jac<1>());

  // SUNMatrix jac = SUNBandMatrix(N_VGetLength(y), 0, 1, sun_ctxt);
  // SUNLinearSolver LS = SUNLinSol_LapackBand(y, jac, sun_ctxt);
  // ARKodeSetLinearSolver(arkode_mem, LS, jac);
  // ARKodeSetJacFn(arkode_mem, create_jac<1>());

  ARKodeSetMaxNonlinIters(arkode_mem, 1000);
  // ARKodeSetNonlinConvCoef(arkode_mem, SUN_RCONST(1.e-4));
  // ARKodeSetSafetyFactor(arkode_mem, 0.8);

  // ARKODE currently approximates the optimal step size with (err)^(-1/(embedded order)), but it's more common to use
  // (err)^(-1/(min(embedded order, order) + 1)). This adjusts the exponent to use the latter.
  ARKodeSetAdaptivityAdjustment(arkode_mem, 0);
  ARKodeSetFixedStepBounds(arkode_mem, 1, 1); // Remove deadzone

  const sunrealtype fac = 1.;
  const sunrealtype reltol = fac * 1.e-2;
  auto abstol = N_VClone(y);
  auto tol_data = N_VGetArrayPointer_Serial(abstol);
  const auto nz = initial_state.t.size();
  for (std::size_t j = 0; j != nz; ++j)
  for (std::size_t j = 0; j != nz; ++j)
  {
    tol_data[j] = fac * 1.e-1;
  }
  for (std::size_t j = 0; j != nz; ++j)
  for (std::size_t j = 0; j != nz; ++j)
  {
    tol_data[nz + j] = fac * 1.e-5;
  }
  for (std::size_t j = 0; j != nz; ++j)
  for (std::size_t j = 0; j != nz; ++j)
  {
    tol_data[2 * nz + j] = fac * 1.e-1;
  }
  for (std::size_t j = 0; j != nz; ++j)
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
      arkode_mem,
      final_time,
      y,
      ARK_NORMAL,
      ARK_ONE_STEP);

  // N_VPrint(y);
  // ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  N_VDestroy(y);
  // SPS: Make RAII wrapper for this.
  ARKodeFree(&arkode_mem);
  SUNMatDestroy(jac);
  SUNLinSolFree(LS);
  return solution;
}
