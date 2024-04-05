#include "sundials_integrator.hpp"
#include "sedimentation.hpp"
#include <cstddef>
#include <vector>
#include <cmath>
#include <iostream>
#include "arkode/arkode_erkstep.h"

SundialsIntegrator::SundialsIntegrator(const RainshaftConstants* constants,
                                       const RainshaftGrid* grid,
                                       const RainshaftProcess* process_exp,
                                       const RainshaftProcess* process_imp,
                                       const RainshaftProcess* process_all,
                                       sundials::Context *sun_ctxt)
  : sun_ctxt(sun_ctxt), user_data{constants, grid, process_exp, process_imp, process_all} {
}

RainshaftDerivedVars SundialsIntegrator::calc_dvars(const RainshaftState& state) const {
  return RainshaftDerivedVars(*user_data.constants, *user_data.grid, state);
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


void tend_jac_to_sun_matrix(const RainshaftTendencyJac& tend_jac, SUNMatrix J) {
  sunindextype nz = SUNDenseMatrix_Rows(J); //length calculation);
  realtype* Jdata = SUNDenseMatrix_Data(J);
  for (sunindextype i = 0; i != nz; ++i) {
    for (sunindextype j = 0; j != nz; ++j) {
      Jdata[j*nz + i] = tend_jac.t_tend_jac[j*nz + i];
    }
  }

  for (sunindextype i = 0; i != nz; ++i) {
    for (sunindextype j = 0; j != nz; ++j) {
      Jdata[j*nz + i] += tend_jac.q_tend_jac[j*nz + i];
    }
  }

  for (sunindextype i = 0; i != nz; ++i) {
    for (sunindextype j = 0; j != nz; ++j) {
      Jdata[j*nz + i] += tend_jac.nr_tend_jac[j*nz + i];
    }
  }

  for (sunindextype i = 0; i != nz; ++i) {
    for (sunindextype j = 0; j != nz; ++j) {
      Jdata[j*nz + i] += tend_jac.qr_tend_jac[j*nz + i];
    }
  }

  // for (sunindextype i = 0; i != nz; ++i) {
  //   for (sunindextype j = 0; j != nz; ++j) {
  //     std::cout << Jdata[j*nz + i]<< " ";
  //   }
  //   std::cout << std::endl;
  // }
  // std::cout << std::endl;
}


int rainshaft_f_imp(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
  // SPS: Should stop using std::vector to reduce copies and allocations.
  RainshaftState state = n_vector_to_state(y);
  RainshaftUserData *cast_data = (RainshaftUserData*) user_data;
  RainshaftDerivedVars dvars = RainshaftDerivedVars(*cast_data->constants,
                                                    *cast_data->grid,
                                                    state);
  RainshaftTendency tend = cast_data->process_imp->calc_tend(*cast_data->constants,
                                                         *cast_data->grid,
                                                         state, dvars);
  tend_to_n_vector(tend, ydot);
  return 0;
}


int rainshaft_f_exp(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
  // SPS: Should stop using std::vector to reduce copies and allocations.
  RainshaftState state = n_vector_to_state(y);
  RainshaftUserData *cast_data = (RainshaftUserData*) user_data;
  RainshaftDerivedVars dvars = RainshaftDerivedVars(*cast_data->constants,
                                                    *cast_data->grid,
                                                    state);
  RainshaftTendency tend = cast_data->process_exp->calc_tend(*cast_data->constants,
                                                         *cast_data->grid,
                                                         state, dvars);
  tend_to_n_vector(tend, ydot);
  return 0;
}


int rainshaft_f_all(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
  // SPS: Should stop using std::vector to reduce copies and allocations.
  RainshaftState state = n_vector_to_state(y);
  RainshaftUserData *cast_data = (RainshaftUserData*) user_data;
  RainshaftDerivedVars dvars = RainshaftDerivedVars(*cast_data->constants,
                                                    *cast_data->grid,
                                                    state);
  RainshaftTendency tend = cast_data->process_all->calc_tend(*cast_data->constants,
                                                         *cast_data->grid,
                                                         state, dvars);
  tend_to_n_vector(tend, ydot);
  return 0;
}


int rainshaft_Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data, 
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  // SPS: Should stop using std::vector to reduce copies and allocations.
  RainshaftState state = n_vector_to_state(y);
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);

  RainshaftUserData *cast_data = (RainshaftUserData*) user_data;
  RainshaftDerivedVars dvars = RainshaftDerivedVars(*cast_data->constants,
                                                    *cast_data->grid,
                                                    state);
  RainshaftTendencyJac tend_jac = cast_data->process_imp->calc_tend_jac(*cast_data->constants,
                                                         *cast_data->grid,
                                                         state, dvars);

  tend_jac_to_sun_matrix(tend_jac, J);
  return 0;
}
