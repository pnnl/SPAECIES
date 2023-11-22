#include "forward_euler_integrator.hpp"
#include "arkode/arkode_erkstep.h"
#include <iostream>

ForwardEulerIntegrator::ForwardEulerIntegrator(const RainshaftConstants* constants,
                                               const RainshaftGrid* grid,
                                               const RainshaftProcess* process,
                                               sundials::Context *sun_ctxt)
  : SundialsIntegrator(constants, grid, process, sun_ctxt) {
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution ForwardEulerIntegrator::integrate(double initial_time,
                                                    double final_time,
                                                    const RainshaftState& initial_state) const {
  std::vector<RainshaftState> states{initial_state};
  std::vector<RainshaftDerivedVars> dvars{calc_dvars(initial_state)};
  sunindextype nz = initial_state.t.size();
  // SPS: It should not assume 4 variables.
  sunindextype num_variables = nz * 4;
  N_Vector y0 = state_to_n_vector(sun_ctxt, initial_state);
  void* arkode_mem = ERKStepCreate(rainshaft_f, initial_time, y0, *sun_ctxt);
  N_VDestroy(y0);
  // Set step size.
  // SPS: Check return value.
  double step_size = final_time - initial_time;
  ERKStepSetFixedStep(arkode_mem, step_size);
  realtype c(0.), a(0.), b(1.);
  ARKodeButcherTable fe_butcher = ARKodeButcherTable_Create(1, 1, 0, &c, &a, &b, NULL);
  // SPS: And this return value.
  ERKStepSetTable(arkode_mem, fe_butcher);
  // SPS: And this return value.
  ERKStepSetUserData(arkode_mem, (void*) &user_data);
  N_Vector yout = N_VNew_Serial(num_variables, *sun_ctxt);
  realtype tret = 0.;
  // SPS: And this return value (and/or returned time).
  ERKStepEvolve(arkode_mem, final_time, yout, &tret, ARK_ONE_STEP);
  RainshaftState new_state = n_vector_to_state(yout);
  N_VDestroy(yout);
  states.push_back(new_state);
  dvars.push_back(calc_dvars(new_state));
  long int num_rhs_evals = 0;
  // SPS: And this return value.
  ERKStepGetNumRhsEvals(arkode_mem, &num_rhs_evals);

  // get last step size that SUNDIALS chose
  double last_step = 0.0;
  ERKStepGetLastStep(arkode_mem, &last_step);

  // current time
  double tcur = 0.0;
  ERKStepGetCurrentTime(arkode_mem, &tcur);

  

  // TODO: 

  // SPS: Make RAII wrapper for this.
  ERKStepFree(&arkode_mem);
  return RainshaftSolution(states, dvars, num_rhs_evals);
}
