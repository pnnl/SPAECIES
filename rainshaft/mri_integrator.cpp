#include "mri_integrator.hpp"
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_mristep.h"
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


MRIIntegrator::MRIIntegrator(const RainshaftConstants* constants,
                                       const int fast_order,
                                       const int slow_order,
                                       const double dt_in,
                                       const double dt_slowfacin,
                                       const RainshaftGrid* grid,
                                       const RainshaftProcess* process_exp,
                                       const RainshaftProcess* process_imp,
                                       const RainshaftProcess* process_all,
                                       sundials::Context *sun_ctxt)
  : SundialsIntegrator(constants, grid, process_exp, process_imp, process_all, sun_ctxt), 
    dt(dt_in), dt_slowfac(dt_slowfacin), fastOrder(fast_order), slowOrder(slow_order) {
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution MRIIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const RainshaftState& initial_state) const {
  int flag;
  std::vector<RainshaftState> states{initial_state};
  std::vector<RainshaftDerivedVars> dvars{calc_dvars(initial_state)};
  sunindextype nz = initial_state.t.size();
  // SPS: It should not assume 4 variables.
  sunindextype num_variables = nz * 4;
  N_Vector y0 = state_to_n_vector(sun_ctxt, initial_state);

  /* use erk_step to step to t=1 */
  void *arkode_mem = NULL;
  arkode_mem = ERKStepCreate(rainshaft_f_all, initial_time, y0, *sun_ctxt);

  ERKStepSetOrder(arkode_mem, 5);
  ERKStepSetUserData(arkode_mem, (void*) &user_data);
  ERKStepSetMaxNumSteps(arkode_mem, 20000000.);

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
  ERKStepSVtolerances(arkode_mem, reltol, abstol);
  N_VDestroy(y0);

  // step to t=1
  N_Vector yout = N_VNew_Serial(num_variables, *sun_ctxt);
  realtype tret = 0.0;
  ERKStepEvolve(arkode_mem, 1.0, yout, &tret, ARK_NORMAL);
  ERKStepFree(&arkode_mem);


  // // step to 1800
  // ERKStepReInit(arkode_mem, rainshaft_f_all, tret, yout);

  // ARKodeButcherTable B = ARKodeButcherTable_Alloc(1, SUNFALSE);
  // B->q = 1;
  // B->p = 0;
  // B->b[0] = 1.0;
  // ERKStepSetTable(arkode_mem, B);

  // // ERKStepSetOrder(arkode_mem, 5);
  // ERKStepSetFixedStep(arkode_mem, dt);


  // auto start = high_resolution_clock::now();
  // while (tret < final_time) {
  //   if (tret + dt > final_time) {
  //     std::cout << "engaging ark_normal" << std::endl;
  //     ERKStepEvolve(arkode_mem, final_time, yout, &tret, ARK_NORMAL);
  //   } else {
  //     ERKStepEvolve(arkode_mem, final_time, yout, &tret, ARK_ONE_STEP);

  //     double last_step = 0.0;
  //     ERKStepGetLastStep(arkode_mem, &last_step);

  //     std::cout << "tret: " << tret << ", Last time step size: " << last_step
  //       << std::endl;
  //     }
  // }
  // auto stop = high_resolution_clock::now();
  // auto duration = duration_cast<milliseconds>(stop - start);
  // cout << duration.count()/1000.0 << endl;




  ////////////////////////////////////////////////////////////////////////////////////

  /* fast (inner) and slow (outer) ARKODE objects */
  void *inner_arkode_mem = NULL;
  void *outer_arkode_mem = NULL;

  /* MRIStepInnerStepper to wrap the inner (fast) ARKStep object */
  MRIStepInnerStepper stepper = NULL;

  /* create an ARKStep object, setting fast (inner) right-hand side
     functions and the initial condition */
  // inner_arkode_mem = ARKStepCreate(NULL, rainshaft_f_imp, initial_time, y0, *sun_ctxt);
  inner_arkode_mem = ARKStepCreate(rainshaft_f_imp, NULL, tret, yout, *sun_ctxt);

  /* setup ARKStep */
  SUNLinearSolver LS = nullptr;  // linear solver memory structure
  SUNMatrix J = nullptr;
  J = SUNDenseMatrix(num_variables, num_variables, *sun_ctxt);
  LS = SUNLinSol_Dense(y0, J, *sun_ctxt);

  // SPS: And this return value.
  int solverOrder = fastOrder;
  ARKStepSetOrder(inner_arkode_mem, solverOrder);

  // SPS: And this return value.
  ARKStepSetUserData(inner_arkode_mem, (void*) &user_data);

  // // Attach linear solver 
  // ARKStepSetLinearSolver(inner_arkode_mem, LS, J);

  // // Jacobian routine
  // ARKStepSetJacFn(inner_arkode_mem, rainshaft_Jac);

  // SPS: And this return value.
  ARKStepSetMaxNumSteps(inner_arkode_mem, 20000000.);
  // ARKStepSetDeduceImplicitRhs(inner_arkode_mem, true);

  // ARKStepSetJacEvalFrequency(arkode_mem, 10);
  // ARKStepSetMaxNonlinIters(inner_arkode_mem, 80);
  // ARKStepSetNonlinConvCoef(inner_arkode_mem, SUN_RCONST(1.e-6));

  // ARKStepSStolerances(arkode_mem, 1.e-6, 1.e-11);  

  // SPS: And this return value.
  ARKStepSVtolerances(inner_arkode_mem, reltol, abstol);

  ARKStepSetSafetyFactor(inner_arkode_mem, 0.85);

  // fixed step for fast solver
  ARKStepSetFixedStep(inner_arkode_mem, dt);

  /* create MRIStepInnerStepper wrapper for the ARKStep memory block */
  ARKStepCreateMRIStepInnerStepper(inner_arkode_mem, &stepper);

  /* create an MRIStep object, setting the slow (outer) right-hand side
     functions and the initial condition */
  outer_arkode_mem = MRIStepCreate(rainshaft_f_exp, NULL, tret, yout, stepper, *sun_ctxt);

  /* Pass udata to user functions */
  MRIStepSetUserData(outer_arkode_mem, (void*) &user_data);

  /* Set the slow step size */
  ARKStepSetFixedStep(inner_arkode_mem, dt);
  MRIStepSetFixedStep(outer_arkode_mem, dt_slowfac*dt);


  MRIStepSVtolerances(outer_arkode_mem, reltol, abstol);

  // second order method
  // ARKodeButcherTable table = ARKodeButcherTable_LoadERK(ARKODE_HEUN_EULER_2_1_2);
  // MRIStepCoupling c = MRIStepCoupling_MIStoMRI(table, 2, 0);
  // ARKodeButcherTable_Free(table);
  // MRIStepSetCoupling(outer_arkode_mem, c);

  MRIStepSetOrder(outer_arkode_mem, slowOrder);

  N_VDestroy(abstol);

  // SPS: And this return value (and/or returned time).
  int i = 0;
  int erkflag;

  // double start;
  std::vector<double> nlsiters, nlsconvfails;
  auto start = high_resolution_clock::now();
  while (tret < final_time) {

    // ARK_NORMAL for last time step so that sundials interpolates
    // back to final_time
    if (tret + dt_slowfac*dt > final_time) {
      // std::cout << "engaging ark_normal" << std::endl;
      MRIStepEvolve(outer_arkode_mem, final_time, yout, &tret, ARK_NORMAL);
    } else {
      MRIStepEvolve(outer_arkode_mem, final_time, yout, &tret, ARK_ONE_STEP);

      // current time
      double tcur = 0.0;
      MRIStepGetCurrentTime(outer_arkode_mem, &tcur);

      double last_step = 0.0;
      MRIStepGetLastStep(outer_arkode_mem, &last_step);

      RainshaftState new_state = n_vector_to_state(yout);
      states.push_back(new_state);
      dvars.push_back(calc_dvars(new_state));

      long int num_rhs_evals_imp = 0;
      long int num_rhs_evals_exp = 0;
      // SPS: And this return value.
      MRIStepGetNumRhsEvals(outer_arkode_mem, &num_rhs_evals_exp, &num_rhs_evals_imp);

      long int num_rhs_evals_imp2 = 0;
      long int num_rhs_evals_exp2 = 0;
      // SPS: And this return value.
      ARKStepGetNumRhsEvals(inner_arkode_mem, &num_rhs_evals_exp2, &num_rhs_evals_imp2);

      // std::cout << "tret: " << tret << ", Last time step size: " << last_step
      //   << ", RHS evals: " << num_rhs_evals_exp+num_rhs_evals_imp+num_rhs_evals_exp2+num_rhs_evals_imp2 << std::endl;
    }


    i++;
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start);
 
  cout << "Slow time step: " << dt_slowfac << "*dt_fast, " << "Fast time step: " << dt << ", Time elapsed: " << duration.count()/1000.0 << endl;


  RainshaftState new_state = n_vector_to_state(yout);
  N_VDestroy(yout);
  states.push_back(new_state);
  dvars.push_back(calc_dvars(new_state));

  long int num_rhs_evals = 0;
  // long int num_rhs_evals_exp = 0;
  // // SPS: And this return value.
  // MRIStepGetNumRhsEvals(outer_arkode_mem, &num_rhs_evals_exp, &num_rhs_evals);

  return RainshaftSolution(states, dvars, num_rhs_evals);
}
