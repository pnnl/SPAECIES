#ifndef RAINSHAFT_INTEGRATOR_HPP
#define RAINSHAFT_INTEGRATOR_HPP
#include "sundials/sundials_context.h"
#include "nvector/nvector_serial.h"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_state.hpp"
#include "rainshaft_tendency.hpp"
#include "rainshaft_process.hpp"
#include "rainshaft_solution.hpp"

class RainshaftIntegrator {

public:

  // Constructor requires timestep.
  RainshaftIntegrator(sundials::Context *sun_ctxt, double dt_in);

  RainshaftSolution integrate(const RainshaftProcess& process,
                              double initial_time,
                              double final_time,
                              const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const RainshaftState& initial_state,
                              const RainshaftDerivedVars& initial_dvars);

  RainshaftState apply_tend(const RainshaftState& state,
                            const RainshaftTendency& tend);

  // SPS: Change pointer type to prevent use of a sun_ctxt that is destroyed.
  sundials::Context *sun_ctxt;

  const double dt;
};

// SPS: See if you can make use of this type safer with pointer changes...
struct RainshaftUserData {
  const RainshaftConstants* constants;
  const RainshaftGrid* grid;
  const RainshaftProcess* process;
};

N_Vector state_to_n_vector(sundials::Context *sun_ctxt, const RainshaftState& state);

void tend_to_n_vector(const RainshaftTendency& tend, N_Vector ydot);

RainshaftState n_vector_to_state(N_Vector y);

int rainshaft_f(realtype t, N_Vector y, N_Vector ydot, void* user_data);

#endif // RAINSHAFT_INTEGRATOR_HPP
