#include "explicit_integrator.hpp"
#include "arkode/arkode_erkstep.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include "sundials/sundials_context.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"
// #include "sunadaptcontroller/sunadaptcontroller_soderlind.h"

#include "sedimentation.hpp"
#include <algorithm>
#include <cmath>

using namespace std;

extern "C" {
  int CFLStabilityFn(N_Vector y, realtype t, realtype *hstab, void *user_data) {
    RainshaftState state = n_vector_to_state(y);
    RainshaftUserData *cast_data = (RainshaftUserData*) user_data;
    const RainshaftConstants *constants = cast_data->constants;
    const RainshaftGrid *grid = cast_data->grid;
    RainshaftDerivedVars dvars = RainshaftDerivedVars(*constants,
                                                      *grid,
                                                      state);
    Sedimentation sed(*constants, false, false);
    double lambdar_top = constants->pi * constants->rhow * constants->nr_top / constants->qr_top;
    lambdar_top = cbrt(lambdar_top);
    std::vector<double> speeds_top = sed.rain_fall_speeds(*constants, constants->rho_top, lambdar_top);
    double cfl_time = dvars.dz[0] / speeds_top[1];
    for (int il = 0; il != grid->nlev - 1; ++il) {
      std::vector<double> speeds = sed.rain_fall_speeds(*constants, dvars.rho_dry[il], dvars.lambdar[il]);
      cfl_time = std::max(cfl_time, dvars.dz[il+1] / speeds[1]);
    }
    *hstab = 0.99 * cfl_time;
    return 0;
  }
}


ExplicitIntegrator::ExplicitIntegrator(const RainshaftConstants* constants,
                                       const double dt_in,
                                       const RainshaftGrid* grid,
                                       const RainshaftProcess* process,
                                       sundials::Context *sun_ctxt)
  : SundialsIntegrator(constants, grid, process, sun_ctxt), dt(dt_in) {
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution ExplicitIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const RainshaftState& initial_state) const {
  int flag;
  std::vector<RainshaftState> states{initial_state};
  std::vector<RainshaftDerivedVars> dvars{calc_dvars(initial_state)};
  sunindextype nz = initial_state.t.size();
  // SPS: It should not assume 4 variables.
  sunindextype num_variables = nz * 4;
  std::cout << "num vars = " << num_variables << std::endl;
  N_Vector y0 = state_to_n_vector(sun_ctxt, initial_state);

  void* arkode_mem = ERKStepCreate(rainshaft_f, initial_time, y0, *sun_ctxt);

  N_VDestroy(y0);

  // SPS: And this return value.
  int solverOrder = 4;
  ERKStepSetOrder(arkode_mem, solverOrder);

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
  

  // // prevent step behavior
  // ERKStepSetFixedStepBounds(arkode_mem, 1.0, 1.0);

  // // // Steven's adaptivity (PID controller, default in arkode)
  // ERKStepSetAdaptivityMethod(arkode_mem, 0, 1, 1, NULL);

  // I-controller
  // ERKStepSetAdaptivityMethod(arkode_mem, 2, 1, 1, NULL);

  // // new adaptivity controls (PID?)
  // SUNAdaptController controller = SUNAdaptController_Soderlind(*sun_ctxt);
  // ERKStepSetAdaptController(arkode_mem, controller);

  // SUNAdaptController_SetParams_PID(controller, 1.0/18.0, 1.0/9.0, 1.0/18.0);

  // SUNAdaptController_SetErrorBias(controller, 1);
  // ERKStepSetErrorBias(arkode_mem, 1);

  ERKStepSetSafetyFactor(arkode_mem, 0.55);


  // SPS: And this return value (and/or returned time).
  realtype tret = 0.;
  double last_step = 0.0;
  int i = 0;

  while (tret < final_time) {

    if (tret + last_step > final_time) {
      ERKStepEvolve(arkode_mem, final_time, yout, &tret, ARK_NORMAL);
    } else {
      ERKStepEvolve(arkode_mem, final_time, yout, &tret, ARK_ONE_STEP);
    }
    
    
    RainshaftState new_state = n_vector_to_state(yout);
    // N_VDestroy(yout);
    states.push_back(new_state);
    dvars.push_back(calc_dvars(new_state));

    long int num_rhs_evals = 0;
    // SPS: And this return value.
    ERKStepGetNumRhsEvals(arkode_mem, &num_rhs_evals);



    // current time
    double tcur = 0.0;
    ERKStepGetCurrentTime(arkode_mem, &tcur);

    long int num_exp_steps = 0;
    ERKStepGetNumExpSteps(arkode_mem, &num_exp_steps);

    long int num_acc_steps = 0;
    ERKStepGetNumAccSteps(arkode_mem, &num_acc_steps);

    long int step_fails = 0;
    ERKStepGetNumErrTestFails(arkode_mem, &step_fails);
 
    ERKStepGetLastStep(arkode_mem, &last_step);

    // // get rel. error vector
    // N_Vector yerr = N_VNew_Serial(NV_LENGTH_S(yout), *sun_ctxt);
    // ERKStepGetEstLocalErrors(arkode_mem, yerr);

    // N_Vector eweight = N_VNew_Serial(NV_LENGTH_S(yout), *sun_ctxt); // N_VClone()
    // ERKStepGetErrWeights(arkode_mem, eweight);

    // N_Vector yerr_ewt = N_VClone(yerr);
    // N_VProd(yerr, eweight, yerr_ewt);

    // // write to file
    // // char myfilename[]
    // FILE * fp;
    // char filepath[256];
    // snprintf (filepath, sizeof(filepath), "/Users/dong9/Desktop/rainshaft_errorvec/rainshaft_errorvec%d.txt", i);

    // fp = fopen (filepath, "w+");
    // N_VPrintFile(yerr_ewt, fp);
    // fclose(fp);

    std::cout << "tret: " << tret << ", Last time step size: " << last_step << ", Acc-limited steps: " << num_acc_steps << ", Error test fails: " << step_fails << std::endl;
    i++;
  }
  


  RainshaftState new_state = n_vector_to_state(yout);
  N_VDestroy(yout);
  states.push_back(new_state);
  dvars.push_back(calc_dvars(new_state));

  long int num_rhs_evals = 0;
  // SPS: And this return value.
  ERKStepGetNumRhsEvals(arkode_mem, &num_rhs_evals);



  FILE *outdata;
  outdata = fopen("sundials_savedresults.csv", "a+");

  ERKStepPrintAllStats(arkode_mem, outdata, SUN_OUTPUTFORMAT_CSV);
  fclose(outdata);

  // SPS: Make RAII wrapper for this.
  ERKStepFree(&arkode_mem);
  return RainshaftSolution(states, dvars, num_rhs_evals);
}
