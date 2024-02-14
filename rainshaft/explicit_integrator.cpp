#include "explicit_integrator.hpp"
#include "arkode/arkode_erkstep.h"
#include "nvector/nvector_serial.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"

ExplicitIntegrator::ExplicitIntegrator(const RainshaftConstants& constants,
                                       const RainshaftGrid& grid,
                                       const RainshaftProcess* const process)
  : SundialsIntegrator(constants, grid, process) {
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution ExplicitIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const RainshaftState& initial_state) const {
  std::vector<RainshaftState> states{initial_state};
  std::vector<RainshaftDerivedVars> dvars{calc_dvars(initial_state)};
  
  auto y = state_to_n_vector(sun_ctxt, initial_state);
  void* arkode_mem = ERKStepCreate(rainshaft_f, initial_time, y, sun_ctxt);
  ERKStepSetUserData(arkode_mem, (void*) &user_data);
  ERKStepSetMaxNumSteps(arkode_mem, -1); // Set no limit on the number of steps

  // Gustafsson’s explicit controller handles stability limited stepsizes well
  auto controller = SUNAdaptController_ExpGus(sun_ctxt);
  // Use no bias and instead rely on a safety factor. Shampine uses these values in the MATLAB ODE suite
  SUNAdaptController_SetErrorBias(controller, 1);
  ERKStepSetSafetyFactor(arkode_mem, 0.8);
  // ARKODE currently approximates the optimal step size with (err)^(-1/(embedded order)), but it's more common to use
  // (err)^(-1/(min(embedded order, order) + 1)). This adjusts the exponent to use the latter.
  ERKStepSetAdaptivityAdjustment(arkode_mem, 0);
  ERKStepSetFixedStepBounds(arkode_mem, 1, 1); // Remove deadzone


  const sunrealtype fac = 1.;
  const sunrealtype reltol = fac * 1.e-2;
  auto abstol = N_VClone(y);
  auto tol_data = N_VGetArrayPointer_Serial(abstol);
  const auto nz = initial_state.t.size();
  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[j] = fac * 1.e-1;
  }
  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[nz + j] = fac * 1.e-5;
  }
  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[2*nz + j] = fac * 1.e-1;
  }
  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[3*nz + j] = fac * 1.e-8;
  }
  ERKStepSVtolerances(arkode_mem, reltol, abstol);
  N_VDestroy(abstol);
  
  sunrealtype tret = 0.;
  ERKStepEvolve(arkode_mem, final_time, y, &tret, ARK_NORMAL);
  const auto new_state = n_vector_to_state(y);
  states.push_back(new_state);
  dvars.push_back(calc_dvars(new_state));


  long int num_rhs_evals = 0;
  ERKStepGetNumRhsEvals(arkode_mem, &num_rhs_evals);

  N_VDestroy(y);
  SUNAdaptController_Destroy(controller);
  // SPS: Make RAII wrapper for this.
  ERKStepFree(&arkode_mem);
  return RainshaftSolution(states, dvars, num_rhs_evals);
}
