#ifndef SUNDIALS_INTEGRATOR_HPP
#define SUNDIALS_INTEGRATOR_HPP
#include "rainshaft_integrator.hpp"
#include "sundials/sundials_context.h"
#include "nvector/nvector_serial.h"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_state.hpp"
#include "rainshaft_process.hpp"
#include "rainshaft_solution.hpp"
#include "rainshaft_tendency.hpp"

// SPS: See if you can make use of this type safer with pointer changes...
struct RainshaftUserData {
  const RainshaftConstants* constants;
  const RainshaftGrid* grid;
  const RainshaftProcess* process;
};

class SundialsIntegrator : public RainshaftIntegrator {

public:

  // SPS: Need to consider wrapping pointer types to check that none of the
  // pointers become invalid while this object exists.
  SundialsIntegrator(const RainshaftConstants* constants,
                     const RainshaftGrid* grid,
                     const RainshaftProcess* process,
                     sundials::Context *sun_ctxt);

  // SPS: Change pointer type to prevent use of a sun_ctxt that is destroyed.
  sundials::Context *sun_ctxt;

protected:

  // SPS: Utility for calculating derived variables.
  RainshaftDerivedVars calc_dvars(const RainshaftState& state) const;

  const RainshaftUserData user_data;

};

N_Vector state_to_n_vector(sundials::Context *sun_ctxt, const RainshaftState& state);

void tend_to_n_vector(const RainshaftTendency& tend, N_Vector ydot);

RainshaftState n_vector_to_state(N_Vector y);

int rainshaft_f(realtype t, N_Vector y, N_Vector ydot, void* user_data);

#endif // SUNDIALS_INTEGRATOR_HPP
