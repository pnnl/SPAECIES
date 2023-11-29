#include "explicit_integrator.hpp"
#include "arkode/arkode_erkstep.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "sundials/sundials_context.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"

#include "sedimentation.hpp"
#include <algorithm>
#include <cmath>

using namespace std;

ExplicitIntegrator::ExplicitIntegrator(const RainshaftConstants* constants,
                                       const RainshaftGrid* grid,
                                       const RainshaftProcess* process,
                                       sundials::Context *sun_ctxt)
  : SundialsIntegrator(constants, grid, process, sun_ctxt) {
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution ExplicitIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const RainshaftState& initial_state) const {
  std::vector<RainshaftState> states{initial_state};
  std::vector<RainshaftDerivedVars> dvars{calc_dvars(initial_state)};
  sunindextype nz = initial_state.t.size();
  // SPS: It should not assume 4 variables.
  sunindextype num_variables = nz * 4;
  N_Vector y0 = state_to_n_vector(sun_ctxt, initial_state);
  //void* arkode_mem = ARKStepCreate(rainshaft_f, NULL, initial_time, y0, *sun_ctxt);

  void* arkode_mem = ERKStepCreate(rainshaft_f, initial_time, y0, *sun_ctxt);

  

  // ERKStepSetFixedStep(arkode_mem, 3.0);

  N_VDestroy(y0);

  // SPS: And this return value.
  ERKStepSetOrder(arkode_mem, 3);

  // custom 2nd order method from Steven
  // sunrealtype a[] = {0,0,0,0,0.3333333333333333,0,0,0,-2.3333333333333335,3,0,0,3,-4,2,0};
  // sunrealtype b[] = {0.5,0,0,0.5};
  // sunrealtype c[] = {0,0.3333333333333333,0.6666666666666666,1};
  // sunrealtype d[] = {0.9523809523809523,0.2222222222222222,0.043478260869565216,-0.2180814354727398};
  // ARKodeButcherTable table = ARKodeButcherTable_Create(4, 2, 1, c, a, b, d);
  // ERKStepSetTable(arkode_mem, table);

  // SPS: And this return value.
  ERKStepSetUserData(arkode_mem, (void*) &user_data);

  // SPS: And this return value.
  ERKStepSetMaxNumSteps(arkode_mem, 20000000.);

  double fac = 1.;
  realtype reltol = fac * 1.e-4;
  N_Vector abstol = N_VNew_Serial(num_variables, *sun_ctxt);
  realtype *tol_data = N_VGetArrayPointer_Serial(abstol);

  // Sean's original tolerances
  // for (sunindextype j = 0; j != nz; ++j) {
  //   tol_data[j] = fac * 1.e-1;
  // }
  // for (sunindextype j = 0; j != nz; ++j) {
  //   tol_data[nz + j] = fac * 1.e-5;
  // }
  // for (sunindextype j = 0; j != nz; ++j) {
  //   tol_data[2*nz + j] = fac * 1.e-1;
  // }
  // for (sunindextype j = 0; j != nz; ++j) {
  //   tol_data[3*nz + j] = fac * 1.e-8;
  // }

  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[j] = 1.e-10;
  }
  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[nz + j] = 1.e-10;
  }
  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[2*nz + j] = 1.e-10;
  }
  for (sunindextype j = 0; j != nz; ++j) {
    tol_data[3*nz + j] = 1.e-10;
  }


  // SPS: And this return value.
  ERKStepSVtolerances(arkode_mem, reltol, abstol);
  N_VDestroy(abstol);
  N_Vector yout = N_VNew_Serial(num_variables, *sun_ctxt);
  

  // prevent step behavior
  ERKStepSetFixedStepBounds(arkode_mem, 1.0, 1.0);

  // // Steven's adaptivity (PID controller, default in arkode)
  // ERKStepSetAdaptivityMethod(arkode_mem, 0, 1, 1, NULL);

  // I-controller
  // ERKStepSetAdaptivityMethod(arkode_mem, 2, 1, 1, NULL);

  // new adaptivity controls (PID?)
  SUNAdaptController controller = SUNAdaptController_Soderlind(*sun_ctxt);
  ERKStepSetAdaptController(arkode_mem, controller);

  SUNAdaptController_SetParams_PID(controller, 1.0/18.0, 1.0/9.0, 1.0/18.0);

  SUNAdaptController_SetErrorBias(controller, 1);
  ERKStepSetErrorBias(arkode_mem, 1);

  ERKStepSetSafetyFactor(arkode_mem, 0.5);


  // set max growth of the adaptive time step (limit it to a factor of 100)
  ERKStepSetMaxFirstGrowth(arkode_mem, 100.0);

  // ERKStepSetMaxStep(arkode_mem, 2.5);
  // ERKStepSetMaxGrowth(arkode_mem, 100.0);


  // SPS: And this return value (and/or returned time).
  realtype tret = 0.;

  ERKStepEvolve(arkode_mem, final_time, yout, &tret, ARK_NORMAL);

  // current time
  double tcur = 0.0;
  ERKStepGetCurrentTime(arkode_mem, &tcur);

  long int num_exp_steps = 0;
  ERKStepGetNumExpSteps(arkode_mem, &num_exp_steps);

  long int num_acc_steps = 0;
  ERKStepGetNumAccSteps(arkode_mem, &num_acc_steps);

  long int step_fails = 0;
  ERKStepGetNumErrTestFails(arkode_mem, &step_fails);

  // get last step size that SUNDIALS chose
  double last_step = 0.0;
  ERKStepGetLastStep(arkode_mem, &last_step);

  std::cout << "tret: " << tret << ", Last time step size: " << last_step << ", Acc-limited steps: " << num_acc_steps << ", Error test fails: " << step_fails << std::endl;

  RainshaftState new_state = n_vector_to_state(yout);
  N_VDestroy(yout);
  states.push_back(new_state);
  dvars.push_back(calc_dvars(new_state));

  long int num_rhs_evals = 0;
  // SPS: And this return value.
  ERKStepGetNumRhsEvals(arkode_mem, &num_rhs_evals);

  // SPS: Make RAII wrapper for this.
  ERKStepFree(&arkode_mem);
  return RainshaftSolution(states, dvars, num_rhs_evals);
}
