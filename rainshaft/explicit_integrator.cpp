#include "explicit_integrator.hpp"
#include "arkode/arkode_erkstep.h"
#include "nvector/nvector_serial.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"

ExplicitIntegrator::ExplicitIntegrator(const RainshaftConstants &constants,
                                       const RainshaftGrid &grid,
                                       const RainshaftProcess *const process,
                                       const double dt,
                                       const int order,
                                       const int steps_per_output)
    : SundialsIntegrator(constants, grid, {process}, steps_per_output),
      dt(dt), order(order)
{
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution ExplicitIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const RainshaftState &initial_state) const
{
  auto y = state_to_n_vector(sun_ctxt, initial_state);
  void *arkode_mem = ERKStepCreate(create_f<0>(), initial_time, y, sun_ctxt);
  ERKStepSetUserData(arkode_mem, (void *)&user_data);
  ERKStepSetFixedStep(arkode_mem, dt);
  ERKStepSetOrder(arkode_mem, order);
  ERKStepSetMaxNumSteps(arkode_mem, -1); // Set no limit on the number of steps
  ERKStepSetStopTime(arkode_mem, final_time);

  // Gustafsson’s explicit controller handles stability limited stepsizes well
  auto controller = SUNAdaptController_I(sun_ctxt);
  // Use no bias and instead rely on a safety factor. Shampine uses these values in the MATLAB ODE suite
  SUNAdaptController_SetErrorBias(controller, 1);
  ERKStepSetSafetyFactor(arkode_mem, 0.9);
  // ARKODE currently approximates the optimal step size with (err)^(-1/(embedded order)), but it's more common to use
  // (err)^(-1/(min(embedded order, order) + 1)). This adjusts the exponent to use the latter.
  ERKStepSetAdaptivityAdjustment(arkode_mem, 0);
  ERKStepSetFixedStepBounds(arkode_mem, 1, 1); // Remove deadzone
  ERKStepSetAdaptController(arkode_mem, controller);

  const sunrealtype fac = 1;
  const sunrealtype reltol = fac * 1.e-14;
  auto abstol = N_VClone(y);
  auto tol_data = N_VGetArrayPointer_Serial(abstol);
  const auto nz = initial_state.t.size();
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
  ERKStepSVtolerances(arkode_mem, reltol, abstol);
  N_VDestroy(abstol);

  auto solution = evolve(
      ERKStepEvolve,
      [arkode_mem]()
      {
        long int num_rhs_evals = 0;
        ERKStepGetNumRhsEvals(arkode_mem, &num_rhs_evals);
        return num_rhs_evals;
      },
      arkode_mem,
      final_time,
      y,
      ARK_NORMAL,
      ARK_ONE_STEP);

  // ERKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  // N_VPrint(y);
  N_VDestroy(y);
  SUNAdaptController_Destroy(controller);
  // SPS: Make RAII wrapper for this.
  ERKStepFree(&arkode_mem);
  return solution;
}
