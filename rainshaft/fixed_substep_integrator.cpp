#include "fixed_substep_integrator.hpp"
#include <cmath>

FixedSubstepIntegrator::FixedSubstepIntegrator(const RainshaftIntegrator *inner_integrator,
                                               double dt_in)
  : dt(dt_in), integrator(inner_integrator) {
}

RainshaftSolution FixedSubstepIntegrator::integrate(double initial_time,
                                                    double final_time,
                                                    const spaecies::State<double>& initial_state) const {
  double time_interval = final_time - initial_time;
  double approx_num_steps = time_interval / dt;
  double rounded_num_steps = std::round(approx_num_steps);
  // Check if a partial time step is needed.
  std::size_t num_full_steps;
  bool need_partial_step = false;
  double partial_step_size = 0.;
  double t_tol = 1.e-8; // SPS: Should make this a setable parameter?
  if (std::abs(approx_num_steps - std::round(approx_num_steps)) <= t_tol) {
    num_full_steps = rounded_num_steps;
  } else {
    num_full_steps = approx_num_steps;
    need_partial_step = true;
    partial_step_size = time_interval - (num_full_steps * dt);
  }
   // empty state
  RainshaftSolution solution(std::vector<spaecies::State<double>>{initial_state}, 0);
  int num_rhs_evals = 0;
  double current_time = initial_time;
  for (std::size_t it = 0; it != num_full_steps; ++it) {
    double next_time = current_time + dt;
    RainshaftSolution step_solution = integrator->integrate(current_time,
                                                       next_time,
                                                       solution.states.back());
    if (it != 0) {
      solution.pop_back();
    }
    step_solution.move_last_to_other(solution);
    num_rhs_evals += step_solution.num_rhs_evals;
    current_time = next_time;
  }
  if (need_partial_step) {
    double next_time = current_time + partial_step_size;
    RainshaftSolution step_solution = integrator->integrate(current_time,
                                                       next_time,
                                                       solution.states.back());
    if (num_full_steps != 0) {
      solution.pop_back();
    }
    step_solution.move_last_to_other(solution);
    num_rhs_evals += step_solution.num_rhs_evals;
  }
  solution.num_rhs_evals = num_rhs_evals;
  return solution;
}
