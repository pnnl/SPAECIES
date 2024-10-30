#include "operator_splitting_integrator.hpp"
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_splittingstep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_stepper.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_lapackdense.h"

OperatorSplittingIntegrator::OperatorSplittingIntegrator(const RainshaftConstants &constants,
                                       const RainshaftGrid &grid,
                                       const RainshaftProcess *const process_exp,
                                       const RainshaftProcess *const process_imp,
                                       const VarDescList& state_descs,
                                       const VarDescList& tend_descs,
                                       const double dt,
                                       const double dt_exp,
                                       const double dt_imp,
                                       const int order,
                                       const int steps_per_output)
    : SundialsIntegrator(constants, grid, {process_exp, process_imp}, state_descs, tend_descs, steps_per_output),
      dt(dt), dt_exp(dt_exp), dt_imp(dt_imp), order(order)
{
}

RainshaftSolution OperatorSplittingIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const StateConst& initial_state) const
{
  const N_Vector y = view_to_n_vector(sun_ctxt, initial_state);

  void *exp_mem = ERKStepCreate(create_f<0>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(exp_mem, (void *)&user_data);
  ARKodeSetFixedStep(exp_mem, dt_exp);
  ARKodeSetOrder(exp_mem, order);
  ARKodeSetMaxNumSteps(exp_mem, -1); // Set no limit on the number of steps
  ARKodeSetStopTime(exp_mem, final_time);

  SUNAdaptController exp_controller = SUNAdaptController_I(sun_ctxt);
  ARKodeSetAdaptController(exp_mem, exp_controller);
  SUNAdaptController_SetErrorBias(exp_controller, 1);
  ARKodeSetSafetyFactor(exp_mem, 0.9);
  ARKodeSetAdaptivityAdjustment(exp_mem, 0);
  ARKodeSetFixedStepBounds(exp_mem, 1, 1);

  void *imp_mem = ARKStepCreate(nullptr, create_f<1>(), initial_time, y, sun_ctxt);
  ARKodeSetUserData(imp_mem, (void *)&user_data);
  ARKodeSetFixedStep(imp_mem, dt_imp);
  ARKodeSetOrder(imp_mem, order);
  ARKodeSetMaxNumSteps(imp_mem, -1); // Set no limit on the number of steps
  ARKodeSetStopTime(imp_mem, final_time);

  SUNMatrix jac = SUNDenseMatrix(N_VGetLength(y), N_VGetLength(y), sun_ctxt);
  SUNLinearSolver LS = SUNLinSol_LapackDense(y, jac, sun_ctxt);
  ARKodeSetLinearSolver(imp_mem, LS, jac);
  ARKodeSetJacFn(imp_mem, create_jac<1>());
  ARKodeSetMaxNonlinIters(imp_mem, 100);
  ARKodeSetDeduceImplicitRhs(imp_mem, true);

  SUNAdaptController imp_controller = SUNAdaptController_I(sun_ctxt);
  ARKodeSetAdaptController(imp_mem, imp_controller);
  SUNAdaptController_SetErrorBias(imp_controller, 1);
  ARKodeSetSafetyFactor(imp_mem, 0.9);
  ARKodeSetAdaptivityAdjustment(imp_mem, 0);
  ARKodeSetFixedStepBounds(imp_mem, 1, 1);

  const sunrealtype fac = 1.;
  const sunrealtype reltol = fac * 1.e-2;
  const N_Vector abstol = N_VClone(y);
  double * const tol_data = N_VGetArrayPointer(abstol);
  const std::size_t nz = user_data.grid.nlev;
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
  ARKodeSVtolerances(exp_mem, reltol, abstol);
  ARKodeSVtolerances(imp_mem, reltol, abstol);
  N_VDestroy(abstol);


  std::array<SUNStepper, 2> steppers{nullptr, nullptr};
  ARKodeCreateSUNStepper(exp_mem, &steppers[0]);
  ARKodeCreateSUNStepper(imp_mem, &steppers[1]);
  void *splitting_mem = SplittingStepCreate(steppers.data(), steppers.size(), initial_time, y, sun_ctxt);
  ARKodeSetFixedStep(splitting_mem, dt);
  ARKodeSetOrder(splitting_mem, order);
  ARKodeSetMaxNumSteps(splitting_mem, -1); // Set no limit on the number of steps
  ARKodeSetStopTime(splitting_mem, final_time);

  const RainshaftSolution solution = evolve(
      ARKodeEvolve,
      [exp_mem, imp_mem]()
      {
        long int num_exp_rhs_evals = 0;
        ERKStepGetNumRhsEvals(exp_mem, &num_exp_rhs_evals);
        long int tmp = 0;
        long int num_imp_rhs_evals = 0;
        ARKStepGetNumRhsEvals(imp_mem, &tmp, &num_imp_rhs_evals);
        return num_exp_rhs_evals + num_imp_rhs_evals;
      },
      []() {},
      splitting_mem,
      final_time,
      y,
      ARK_NORMAL,
      ARK_ONE_STEP);

  N_VDestroy(y);
  SUNAdaptController_Destroy(exp_controller);
  SUNAdaptController_Destroy(imp_controller);
  SUNStepper_Destroy(&steppers[0]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&exp_mem);
  ARKodeFree(&imp_mem);
  ARKodeFree(&splitting_mem);
  return solution;
}
