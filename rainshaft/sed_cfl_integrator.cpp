#include "sed_cfl_integrator.hpp"
#include "rainshaft_derived_vars.hpp"
#include <algorithm>
#include <cmath>

SedCflIntegrator::SedCflIntegrator(const RainshaftConstants* constants,
                                   const RainshaftGrid* grid,
                                   const VarDescList& tend_descs,
                                   const Sedimentation *sedimentation)
  : constants(constants), grid(grid), sed(sedimentation), tend_descs(tend_descs) {
}

RainshaftSolution SedCflIntegrator::integrate(double initial_time,
                                              double final_time,
                                              const StateConst& initial_state) const {
  int num_rhs_evals = 0;
  double current_time = initial_time;
  double time_remaining = final_time - current_time;
  std::size_t it = 0;
  State temp_state(initial_state.deep_copy());
  double *state_data = temp_state.data();
  Tendency tend(tend_descs);
  double *tend_data = tend.data();
  while (time_remaining > 0) {
    RainshaftDerivedVars dvar = RainshaftDerivedVars(*constants, *grid, temp_state);
    double lambdar_top = std::cbrt(constants->pi * constants->rhow
                                   * constants->nr_top / constants->qr_top);
    const auto [v0, v3] = sed->rain_fall_speeds(*constants, dvar.rho_dry[0],
                                        lambdar_top);
    double max_cfl_num = time_remaining * v3 / dvar.dz[0];
    for (std::size_t il = 0; il != grid->nlev; ++il) {
      const auto [v0, v3] = sed->rain_fall_speeds(*constants, dvar.rho_dry[il],
                                     dvar.lambdar[il]);
      max_cfl_num = std::max(max_cfl_num,
                             time_remaining * v3 / dvar.dz[il]);
    }
    std::size_t estimated_steps_left = ((std::size_t) max_cfl_num) + 1;
    double step_size = time_remaining / estimated_steps_left;
    double next_time = std::min(current_time + step_size, final_time);
    // Take a forward Euler step.
    sed->calc_tend(*constants, *grid, temp_state, dvar, tend);
    for (std::size_t i = 0; i != tend.size(); ++i) {
      state_data[i] = std::max(state_data[i] + step_size * tend_data[i], 0.);
    }
    ++num_rhs_evals;
    current_time = next_time;
    time_remaining = final_time - current_time;
    ++it;
  }
  return RainshaftSolution({temp_state}, num_rhs_evals);
}
