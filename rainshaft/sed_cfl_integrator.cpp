#include "sed_cfl_integrator.hpp"
#include "rainshaft_derived_vars.hpp"
#include <cmath>

SedCflIntegrator::SedCflIntegrator(const RainshaftConstants* constants,
                                   const RainshaftGrid* grid,
                                   const Sedimentation *sedimentation,
                                   const RainshaftIntegrator *integrator)
  : constants(constants), grid(grid), sed(sedimentation), sed_int(integrator) {
}

RainshaftSolution SedCflIntegrator::integrate(double initial_time,
                                              double final_time,
                                              const StateConst& initial_state) const {
  RainshaftSolution solution({initial_state}, 0);
  int num_rhs_evals = 0;
  double current_time = initial_time;
  double time_remaining = final_time - current_time;
  std::size_t it = 0;
  while (time_remaining > 0) {
    RainshaftDerivedVars dvar = RainshaftDerivedVars(*constants, *grid, solution.states.back());
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
    RainshaftSolution step_solution = sed_int->integrate(current_time,
                                                    next_time,
                                                    solution.states.back());
    num_rhs_evals += step_solution.num_rhs_evals;
    solution.pop_back();
    step_solution.move_last_to_other(solution);
    current_time = next_time;
    time_remaining = final_time - current_time;
    ++it;
  }
  solution.num_rhs_evals = num_rhs_evals;
  return solution;
}
