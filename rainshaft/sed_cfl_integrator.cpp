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
                                              const RainshaftState& initial_state) const {
  std::vector<RainshaftState> states{initial_state};
  std::vector<RainshaftDerivedVars> dvars{RainshaftDerivedVars(*constants, *grid, initial_state)};
  int num_rhs_evals = 0;
  double current_time = initial_time;
  double time_remaining = final_time - current_time;
  std::size_t it = 0;
  while (time_remaining > 0) {
    double lambdar_top = std::cbrt(constants->pi * constants->rhow
                                   * constants->nr_top / constants->qr_top);
    auto speeds = sed->rain_fall_speeds(*constants, dvars[it].rho_dry[0],
                                        lambdar_top);
    double v3 = speeds[1];
    double max_cfl_num = time_remaining * v3 / dvars[it].dz[0];
    for (std::size_t il = 0; il != grid->nlev; ++il) {
      speeds = sed->rain_fall_speeds(*constants, dvars[it].rho_dry[il],
                                     dvars[it].lambdar[il]);
      v3 = speeds[1];
      max_cfl_num = std::max(max_cfl_num,
                             time_remaining * v3 / dvars[it].dz[il]);
    }
    std::size_t estimated_steps_left = ((std::size_t) max_cfl_num) + 1;
    double step_size = time_remaining / estimated_steps_left;
    double next_time = std::min(current_time + step_size, final_time);
    RainshaftSolution solution = sed_int->integrate(current_time,
                                                    next_time,
                                                    states[it]);
    states.push_back(solution.states.back());
    dvars.push_back(solution.dvars.back());
    num_rhs_evals += solution.num_rhs_evals;
    current_time = next_time;
    time_remaining = final_time - current_time;
    ++it;
  }
  return RainshaftSolution(states, dvars, num_rhs_evals);
}
