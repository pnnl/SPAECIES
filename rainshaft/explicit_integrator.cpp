#include "explicit_integrator.hpp"
#include "arkode/arkode_erkstep.h"

ExplicitIntegrator::ExplicitIntegrator(sundials::Context *sun_ctxt)
  : RainshaftIntegrator(sun_ctxt) {
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution ExplicitIntegrator::integrate(const RainshaftProcess& process,
                                                 double initial_time,
                                                 double final_time,
                                                 const RainshaftConstants& constants,
                                                 const RainshaftGrid& grid,
                                                 const RainshaftState& initial_state,
                                                 const RainshaftDerivedVars& initial_dvars) {
  std::vector<RainshaftState> states{initial_state};
  std::vector<RainshaftDerivedVars> dvars{initial_dvars};
  sunindextype nz = initial_state.t.size();
  // SPS: It should not assume 4 variables.
  sunindextype num_variables = nz * 4;
  // SPS: initial_dvars unused, find a way to use it or remove it as an argument.
  N_Vector y0 = state_to_n_vector(sun_ctxt, initial_state);
  // Make userdata object.
  RainshaftUserData user_data = {&constants, &grid, &process};
  //void* arkode_mem = ARKStepCreate(rainshaft_f, NULL, initial_time, y0, *sun_ctxt);
  void* arkode_mem = ERKStepCreate(rainshaft_f, initial_time, y0, *sun_ctxt);
  // SPS: And this return value.
  ERKStepSetOrder(arkode_mem, 2);
  // SPS: And this return value.
  ERKStepSetUserData(arkode_mem, (void*) &user_data);
  // SPS: And this return value.
  ERKStepSetMaxNumSteps(arkode_mem, 20000000.);
  realtype reltol = 1.e-2;
  N_Vector abstol = N_VNew_Serial(num_variables, *sun_ctxt);
  realtype *tol_data = N_VGetArrayPointer_Serial(abstol);
  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[j] = 1.e-1;
  }
  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[nz + j] = 1.e-6;
  }
  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[2*nz + j] = 1.e-1;
  }
  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[3*nz + j] = 1.e-8;
  }
  // SPS: And this return value.
  ERKStepSVtolerances(arkode_mem, reltol, abstol);
  N_Vector yout = N_VNew_Serial(num_variables, *sun_ctxt);
  realtype tret = 0.;
  // SPS: And this return value (and/or returned time).
  ERKStepEvolve(arkode_mem, final_time, yout, &tret, ARK_NORMAL);
  RainshaftState new_state = n_vector_to_state(yout);
  states.push_back(new_state);
  dvars.push_back(RainshaftDerivedVars(constants, grid, new_state));
  long int num_rhs_evals = 0, num_implicit_evals = 0;
  // SPS: And this return value.
  //ARKStepGetNumRhsEvals(arkode_mem, &num_rhs_evals,
  //                      &num_implicit_evals);
  //num_rhs_evals += num_implicit_evals;
  ERKStepGetNumRhsEvals(arkode_mem, &num_rhs_evals);

  // SPS: Make RAII wrapper for this.
  ERKStepFree(&arkode_mem);
  return RainshaftSolution(states, dvars, num_rhs_evals);
}
