#include "sed_cfl_integrator.hpp"
#include "rainshaft_derived_vars.hpp"
#include <algorithm>
#include <cmath>

SedCflIntegrator::SedCflIntegrator(const RainshaftConstants& constants,
                                   const RainshaftGrid& grid,
                                   const SizeLimiters& size_limiters,
                                   const VarDescList& tend_descs,
                                   const Sedimentation& sedimentation,
                                   const bool regularize_lambdar)
  : constants(constants), grid(grid), size_limiters(size_limiters), sed(sedimentation), tend_descs(tend_descs),
    regularize_lambdar(regularize_lambdar) {
}

RainshaftSolution SedCflIntegrator::integrate(double initial_time,
                                              double final_time,
                                              const StateConst& initial_state,
                                              int& error_flag) const {
  int num_rhs_evals = 0;
  double current_time = initial_time;
  double time_remaining = final_time - current_time;
  std::size_t it = 0;
  State temp_state(initial_state.deep_copy());
  double *state_data = temp_state.data();
  Tendency tend(tend_descs);
  double *tend_data = tend.data();
  VarMut nr = temp_state.get_variable("nr");
  VarConst qr = temp_state.get_variable("qr");
  while (time_remaining > 0) {
    std::fill_n(tend_data, tend.size(), 0.0);
    RainshaftDerivedVars dvars = RainshaftDerivedVars(constants, grid, temp_state, regularize_lambdar);
    const double max_cfl_num = time_remaining / sed.calc_max_step(constants, grid, dvars);
    const std::size_t estimated_steps_left = ((std::size_t) max_cfl_num) + 1;
    const double step_size = time_remaining / estimated_steps_left;
    // Take a forward Euler step.
    sed.calc_tend(constants, grid, temp_state, dvars, tend);
    // Apply tendency, making sure all values stay above 0.
    for (std::size_t i = 0; i != tend.size(); ++i) {
      state_data[i] = std::max(state_data[i] + step_size * tend_data[i], 0.);
    }
    // Apply size limiters.
    for (std::size_t i = 0; i != grid.nlev; ++i) {
      nr[i] = size_limiters.limited_nr(nr[i], qr[i]);
    }
    ++num_rhs_evals;
    current_time = std::min(current_time + step_size, final_time);
    time_remaining = final_time - current_time;
    ++it;
  }
  return RainshaftSolution({temp_state}, num_rhs_evals);
}
