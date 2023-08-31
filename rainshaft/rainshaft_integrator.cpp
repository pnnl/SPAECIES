#include "rainshaft_integrator.hpp"
#include <cstddef>
#include <vector>
#include "arkode/arkode_erkstep.h"

RainshaftIntegrator::RainshaftIntegrator(sundials::Context *sun_ctxt)
  : sun_ctxt(sun_ctxt) {
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
