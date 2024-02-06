#include "sundials_integrator.hpp"
#include "nvector/nvector_serial.h"
#include <stdexcept>
#include <sstream>
#include <cstddef>
#include <vector>

SundialsIntegrator::SundialsIntegrator(const RainshaftConstants& constants,
                                       const RainshaftGrid& grid,
                                       const RainshaftProcess* const process)
  : user_data{constants, grid, process} {
    SUNContext_PushErrHandler(sun_ctxt, [](int line, const char *func, const char *file, const char *msg,
      SUNErrCode err_code, void *err_user_data, SUNContext sunctx)
    {
      std::ostringstream buf;
      buf << msg << "\nIn " << func << " on line " << line << " of " << file;
      throw std::logic_error(buf.str());
    }, nullptr);
}

RainshaftDerivedVars SundialsIntegrator::calc_dvars(const RainshaftState& state) const {
  return RainshaftDerivedVars(user_data.constants, user_data.grid, state);
}

N_Vector state_to_n_vector(const sundials::Context& sun_ctxt, const RainshaftState& state) {
  sunindextype nz = state.t.size();
  // SPS: It should not assume 4 variables; move some of this to RainshaftState itself?
  sunindextype num_variables = nz * 4;
  N_Vector y = N_VNew_Serial(num_variables, sun_ctxt);
  sunrealtype *ydata = N_VGetArrayPointer_Serial(y);
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
  sunrealtype *ydata = N_VGetArrayPointer_Serial(ydot);
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
  sunrealtype* ydata = N_VGetArrayPointer(y);
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

int rainshaft_f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data) {
  // SPS: Should stop using std::vector to reduce copies and allocations.
  RainshaftState state = n_vector_to_state(y);
  RainshaftUserData *cast_data = (RainshaftUserData*) user_data;
  RainshaftDerivedVars dvars = RainshaftDerivedVars(cast_data->constants,
                                                    cast_data->grid,
                                                    state);
  RainshaftTendency tend = cast_data->process->calc_tend(cast_data->constants,
                                                         cast_data->grid,
                                                         state, dvars);
  tend_to_n_vector(tend, ydot);
  return 0;
}
