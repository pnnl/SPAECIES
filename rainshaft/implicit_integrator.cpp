#include "implicit_integrator.hpp"
#include "arkode/arkode_arkstep.h"
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
using namespace std::chrono;


ImplicitIntegrator::ImplicitIntegrator(const RainshaftConstants* constants,
                                       const double dt_in,
                                       const RainshaftGrid* grid,
                                       const RainshaftProcess* process_exp,
                                       const RainshaftProcess* process_imp,
                                       sundials::Context *sun_ctxt)
  : SundialsIntegrator(constants, grid, process_exp, process_imp, sun_ctxt), dt(dt_in) {
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution ImplicitIntegrator::integrate(double initial_time,
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

  // void* arkode_mem = ARKStepCreate(rainshaft_f_imp, NULL, initial_time, y0, *sun_ctxt);
  void* arkode_mem = ARKStepCreate(NULL, rainshaft_f_imp, initial_time, y0, *sun_ctxt);

  SUNLinearSolver LS = nullptr;  // linear solver memory structure
  SUNMatrix J = nullptr;
  J = SUNDenseMatrix(num_variables, num_variables, *sun_ctxt);
  LS = SUNLinSol_Dense(y0, J, *sun_ctxt);

  N_VDestroy(y0);

  // SPS: And this return value.
  int solverOrder = 5;
  ARKStepSetOrder(arkode_mem, solverOrder);

  // SPS: And this return value.
  ARKStepSetUserData(arkode_mem, (void*) &user_data);

  // Attach linear solver 
  ARKStepSetLinearSolver(arkode_mem, LS, J);

  // // Jacobian routine
  // ARKStepSetJacFn(arkode_mem, rainshaft_Jac);

  // SPS: And this return value.
  ARKStepSetMaxNumSteps(arkode_mem, 20000000.);

  ARKStepSetMaxNonlinIters(arkode_mem, 80);
  ARKStepSetNonlinConvCoef(arkode_mem, SUN_RCONST(1.e-6));

  // ARKStepSStolerances(arkode_mem, 1.e-6, 1.e-11);  

  double fac = 1.;
  realtype reltol = fac * 1.e-4;
  N_Vector abstol = N_VNew_Serial(num_variables, *sun_ctxt);
  realtype *tol_data = N_VGetArrayPointer_Serial(abstol);

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
  ARKStepSVtolerances(arkode_mem, reltol, abstol);
  N_VDestroy(abstol);
  N_Vector yout = N_VNew_Serial(num_variables, *sun_ctxt);

  ARKStepSetSafetyFactor(arkode_mem, 0.55);


  // SPS: And this return value (and/or returned time).
  realtype tret = 0.;
  int i = 0;

  // double start;
  auto start = high_resolution_clock::now();
  while (tret < final_time) {

    // use adaptive implicit stepper until t = 312.0
    // ReInit here is in case of switching to a different method
    if (tret > 317.507 & tret < 317.507 + dt) { // probably a better way to do this range
                                       // so that we don't call ReInit multiple times
      ARKStepReInit(arkode_mem, NULL, rainshaft_f_imp, tret, yout);

      ARKStepSetOrder(arkode_mem, 5);


      ARKStepSetJacFn(arkode_mem, rainshaft_Jac);

      // ARKStepSetMaxNonlinIters(arkode_mem, 160);
      // ARKStepSetNonlinConvCoef(arkode_mem, SUN_RCONST(1.e-10));
      std::cout << "Switching method" << std::endl;
      start = high_resolution_clock::now();
    }

    // use fixed step for t > 312
    if (tret > 312.0) {
      ARKStepSetFixedStep(arkode_mem, dt);
    } 


    // ARK_NORMAL for last time step so that sundials interpolates
    // back to final_time
    if (tret + 100.0 > final_time) {
      ARKStepEvolve(arkode_mem, final_time, yout, &tret, ARK_NORMAL);
    } else {
      ARKStepEvolve(arkode_mem, final_time, yout, &tret, ARK_ONE_STEP);
    }
    
    
    RainshaftState new_state = n_vector_to_state(yout);
    // N_VDestroy(yout);
    states.push_back(new_state);
    dvars.push_back(calc_dvars(new_state));

    long int num_rhs_evals = 0;
    long int num_rhs_evals_exp = 0;
    // SPS: And this return value.
    ARKStepGetNumRhsEvals(arkode_mem, &num_rhs_evals_exp, &num_rhs_evals);



    // current time
    double tcur = 0.0;
    ARKStepGetCurrentTime(arkode_mem, &tcur);

    long int num_exp_steps = 0;
    ARKStepGetNumExpSteps(arkode_mem, &num_exp_steps);

    long int num_acc_steps = 0;
    ARKStepGetNumAccSteps(arkode_mem, &num_acc_steps);

    long int step_fails = 0;
    ARKStepGetNumErrTestFails(arkode_mem, &step_fails);

    double last_step = 0.0;
    ARKStepGetLastStep(arkode_mem, &last_step);

    long int nniters = 0;
    ARKStepGetNumNonlinSolvIters(arkode_mem, &nniters);

    long int nncfails = 0;
    ARKStepGetNumNonlinSolvConvFails(arkode_mem, &nncfails);

    // get rel error vector
    N_Vector yerr = N_VNew_Serial(NV_LENGTH_S(yout), *sun_ctxt);
    ARKStepGetEstLocalErrors(arkode_mem, yerr);

    N_Vector eweight = N_VNew_Serial(NV_LENGTH_S(yout), *sun_ctxt); // N_VClone()
    ARKStepGetErrWeights(arkode_mem, eweight);

    N_Vector yerr_ewt = N_VClone(yerr);
    N_VProd(yerr, eweight, yerr_ewt);

    // // write to file
    // // char myfilename[]
    // FILE * fp;
    // char filepath[256];
    // snprintf (filepath, sizeof(filepath), "/Users/dong9/Desktop/rainshaft_errorvec/rainshaft_errorvec%d.txt", i);

    // fp = fopen (filepath, "w+");
    // N_VPrintFile(yerr_ewt, fp);
    // // myfile.close();
    // fclose(fp);


    // // write to file
    // // char myfilename[]
    // ARKStepGetJac(arkode_mem, &J);

    // char filepath2[256];
    // snprintf (filepath2, sizeof(filepath2), "/Users/dong9/Desktop/Jacobians/jacobian%d.txt", i);

    // fp = fopen (filepath2, "w+");
    // SUNDenseMatrix_Print(J, fp);
    // fclose(fp);

    // realtype* Jdata = SUNDenseMatrix_Data(J);
    // for (sunindextype i = num_variables/2; i != num_variables; ++i) {
    //   for (sunindextype j = num_variables/2; j != num_variables; ++j) {
    //     std::cout << Jdata[j*num_variables + i]<< " ";
    //   }
    //   std::cout << std::endl;
    // }
    // std::cout << std::endl;

    std::cout << "tret: " << tret << ", Last time step size: " << last_step << ", Acc-limited steps: " << num_acc_steps << ", Error test fails: " << step_fails << ", NLS iters: " << nniters << ", NLS conv. fails: " << nncfails << std::endl;
    i++;
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start);
 
  // To get the value of duration use the count()
  // member function on the duration object
  cout << duration.count()/1000.0 << endl;


  RainshaftState new_state = n_vector_to_state(yout);
  N_VDestroy(yout);
  states.push_back(new_state);
  dvars.push_back(calc_dvars(new_state));

  long int num_rhs_evals = 0;
  long int num_rhs_evals_exp = 0;
  // SPS: And this return value.
  ARKStepGetNumRhsEvals(arkode_mem, &num_rhs_evals_exp, &num_rhs_evals);



  FILE *outdata;
  outdata = fopen("sundials_savedresults.csv", "a+");

  ARKStepPrintAllStats(arkode_mem, outdata, SUN_OUTPUTFORMAT_CSV);
  fclose(outdata);

  // SPS: Make RAII wrapper for this.
  ARKStepFree(&arkode_mem);
  return RainshaftSolution(states, dvars, num_rhs_evals);
}
