#include "rainshaft_integrator.hpp"
#include <cstddef>
#include <iostream>
#include <vector>
#include "arkode/arkode_erkstep.h"

RainshaftIntegrator::RainshaftIntegrator(sundials::Context *sun_ctxt, double dt_in)
  : sun_ctxt(sun_ctxt), dt(dt_in) {
}

// SPS: Need to generalize this to get output states at arbitary times.
RainshaftSolution RainshaftIntegrator::integrate(const RainshaftProcess& process,
                                                 double initial_time,
                                                 double final_time,
                                                 const RainshaftConstants& constants,
                                                 const RainshaftGrid& grid,
                                                 const RainshaftState& initial_state,
                                                 const RainshaftDerivedVars& initial_dvars) {
  std::size_t num_steps = (final_time - initial_time) / dt;
  std::vector<RainshaftState> states{initial_state};
  std::vector<RainshaftDerivedVars> dvars{initial_dvars};
  sunindextype nz = initial_state.t.size();
  // SPS: It should not assume 4 variables.
  sunindextype num_variables = nz * 4;
  // SPS: initial_dvars unused, find a way to use it or remove it as an argument.
  N_Vector y0 = state_to_n_vector(sun_ctxt, initial_state);
  // Make userdata object.
  RainshaftUserData user_data = {&constants, &grid, &process};
  void* arkode_mem = ERKStepCreate(rainshaft_f, initial_time, y0, *sun_ctxt);
  // SPS: Check return value.
  ERKStepSetFixedStep(arkode_mem, dt);
  realtype c(0.), a(0.), b(1.);
  ARKodeButcherTable butcher = ARKodeButcherTable_Create(1, 1, 0, &c, &a, &b, NULL);
  // SPS: And this return value.
  ERKStepSetTable(arkode_mem, butcher);
  // SPS: And this return value.
  ERKStepSetUserData(arkode_mem, (void*) &user_data);
  // SPS: And this return value.
  ERKStepSetMaxNumSteps(arkode_mem, num_steps + 1);
  N_Vector yout = N_VNew_Serial(num_variables, *sun_ctxt);
  realtype tret = 0.;
  for (std::size_t it = 1; it != num_steps + 1; ++it) {
    // SPS: And this return value (and/or returned time).
    ERKStepEvolve(arkode_mem, dt * it, yout, &tret, ARK_ONE_STEP);
    RainshaftState new_state = n_vector_to_state(yout);
    states.push_back(new_state);
    dvars.push_back(RainshaftDerivedVars(constants, grid, new_state));
  }
  //for (std::size_t it = 0; it != num_steps; ++it) {
  //  RainshaftTendency tend = process.calc_tend(constants, grid, states[it], dvars[it]);
  //  RainshaftState new_state = apply_tend(states[it], tend);
  //}

  // SPS: Make RAII wrapper for this.
  ERKStepFree(&arkode_mem);
  return RainshaftSolution(states, dvars);
}

N_Vector state_to_n_vector(sundials::Context *sun_ctxt, const RainshaftState& state) {
  sunindextype nz = state.t.size();
  // SPS: It should not assume 4 variables; move some of this to RainshaftState itself?
  sunindextype num_variables = nz * 4;
  N_Vector y = N_VNew_Serial(num_variables, *sun_ctxt);
  realtype *ydata = N_VGetArrayPointer_Serial(y);
  for (sunindextype j = 0; j != nz; ++j) {
    ydata[j] = state.t[j];
  }
  for (sunindextype j = 0; j != nz; ++j) {
    ydata[nz + j] = state.q[j];
  }
  for (sunindextype j = 0; j != nz; ++j) {
    ydata[2*nz + j] = state.nr[j];
  }
  for (sunindextype j = 0; j != nz; ++j) {
    ydata[3*nz + j] = state.qr[j];
  }
  return y;
}

// SPS: should reuse most of the equivalent state function here
void tend_to_n_vector(const RainshaftTendency& tend, N_Vector ydot) {
  sunindextype nz = tend.t_tend.size();
  // SPS: still bad to assume these 4 variables.
  sunindextype num_variables = nz * 4;
  realtype *ydata = N_VGetArrayPointer_Serial(ydot);
  for (sunindextype j = 0; j != nz; ++j) {
    ydata[j] = tend.t_tend[j];
  }
  for (sunindextype j = 0; j != nz; ++j) {
    ydata[nz + j] = tend.q_tend[j];
  }
  for (sunindextype j = 0; j != nz; ++j) {
    ydata[2*nz + j] = tend.nr_tend[j];
  }
  for (sunindextype j = 0; j != nz; ++j) {
    ydata[3*nz + j] = tend.qr_tend[j];
  }
}

RainshaftState n_vector_to_state(N_Vector y) {
  sunindextype nz = N_VGetLength(y) / 4;
  realtype* ydata = N_VGetArrayPointer(y);
  std::vector<double> t(nz), q(nz), nr(nz), qr(nz);
  for (sunindextype j = 0; j != nz; ++j) {
    t[j] = ydata[j];
  }
  for (sunindextype j = 0; j != nz; ++j) {
    q[j] = ydata[nz + j];
  }
  for (sunindextype j = 0; j != nz; ++j) {
    nr[j] = ydata[2*nz + j];
  }
  for (sunindextype j = 0; j != nz; ++j) {
    qr[j] = ydata[3*nz + j];
  }
  return RainshaftState(t, q, nr, qr);
}

int rainshaft_f(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
  // SPS: Should stop using std::vector to reduce copies and allocations.
  RainshaftState state = n_vector_to_state(y);
  RainshaftUserData *cast_data = (RainshaftUserData*) user_data;
  RainshaftDerivedVars dvars = RainshaftDerivedVars(*cast_data->constants,
                                                    *cast_data->grid,
                                                    state);
  RainshaftTendency tend = cast_data->process->calc_tend(*cast_data->constants,
                                                         *cast_data->grid,
                                                         state, dvars);
  tend_to_n_vector(tend, ydot);
  return 0;
}

RainshaftState RainshaftIntegrator::apply_tend(const RainshaftState& state,
                                               const RainshaftTendency& tend) {
  std::vector<double> t, q, nr, qr;
  for (std::size_t i = 0; i != state.t.size(); ++i) {
    t.push_back(state.t[i] + dt*tend.t_tend[i]);
    q.push_back(state.q[i] + dt*tend.q_tend[i]);
    nr.push_back(state.nr[i] + dt*tend.nr_tend[i]);
    qr.push_back(state.qr[i] + dt*tend.qr_tend[i]);
  }
  return RainshaftState(t, q, nr, qr);
}
