#ifndef SUNDIALS_INTEGRATOR_HPP
#define SUNDIALS_INTEGRATOR_HPP
#include "rainshaft_integrator.hpp"
#include "sundials/sundials_context.hpp"
#include "sundials/sundials_nvector.h"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_state.hpp"
#include "rainshaft_process.hpp"
#include "rainshaft_solution.hpp"
#include "rainshaft_tendency.hpp"

struct RainshaftUserData {
  const RainshaftConstants& constants;
  const RainshaftGrid& grid;
  const RainshaftProcess* const process;
};

class SundialsIntegrator : public RainshaftIntegrator {

public:
  SundialsIntegrator(const RainshaftConstants& constants,
                     const RainshaftGrid& grid,
                     const RainshaftProcess* const process);

  const sundials::Context sun_ctxt;

protected:

  // SPS: Utility for calculating derived variables.
  RainshaftDerivedVars calc_dvars(const RainshaftState& state) const;

  const RainshaftUserData user_data;

};

N_Vector state_to_n_vector(const sundials::Context& sun_ctxt, const RainshaftState& state);

void tend_to_n_vector(const RainshaftTendency& tend, N_Vector ydot);

RainshaftState n_vector_to_state(N_Vector y);

int rainshaft_f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

#endif // SUNDIALS_INTEGRATOR_HPP
