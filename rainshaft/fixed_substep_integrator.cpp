#include "fixed_substep_integrator.hpp"
#include <cmath>
#include <iostream>

FixedSubstepIntegrator::FixedSubstepIntegrator(const RainshaftIntegrator *inner_integrator,
                                               double dt_in)
  : dt(dt_in), integrator(inner_integrator) {
}

RainshaftSolution FixedSubstepIntegrator::integrate(double initial_time,
                                                    double final_time,
                                                    const RainshaftState& initial_state) const {
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
  std::vector<RainshaftState> states{initial_state};
  std::vector<RainshaftDerivedVars> dvars;
  int num_rhs_evals = 0;
  double current_time = initial_time;


  double next_time = current_time + dt;
  RainshaftSolution solution = integrator->integrate(current_time,
                                                     next_time,
                                                     states.back());

  for (std::size_t it = 0; it < solution.states.size(); ++it) {
    dvars.push_back(solution.dvars[it]);

    if (it != 0) {
      states.push_back(solution.states[it]);
    }

    double next_time = current_time + dt;
    num_rhs_evals += solution.num_rhs_evals;
    current_time = next_time;
  }


  // for (std::size_t it = 0; it != num_full_steps; ++it) {
  //   double next_time = current_time + dt;
  //   RainshaftSolution solution = integrator->integrate(current_time,
  //                                                      next_time,
  //                                                      states.back());
  //   // This class currently doesn't store constants/grid, so we rely
  //   // on the integrator down a level to calculate the initial dvars
  //   // for us.
  //   // SPS: Does this need fixed? Do we need to handle the case of
  //   // zero total steps, for instance?

  //   // if (it == 0) {
  //   //   dvars.push_back(solution.dvars[0]);
  //   // } else {
  //   //   states.pop_back();
  //   // }
  //   // states.push_back(solution.states.back());
  //   // if ((it == num_full_steps - 1) && !need_partial_step) {
  //   //   dvars.push_back(solution.dvars.back());
  //   // }


  //   dvars.push_back(solution.dvars[0]);

  //   if (it != 0) {
  //     states.push_back(solution.states.back());
  //   }

  //   num_rhs_evals += solution.num_rhs_evals;
  //   current_time = next_time;
  // }


  if (need_partial_step) {
    double next_time = current_time + partial_step_size;
    RainshaftSolution solution = integrator->integrate(current_time,
                                                       next_time,
                                                       states.back());
    // Push the initial dvars in the case that the partial step is the
    // *only* step.
    if (num_full_steps == 0) {
      dvars.push_back(solution.dvars[0]);
    } else {
      states.pop_back();
    }
    states.push_back(solution.states.back());
    dvars.push_back(solution.dvars.back());
    num_rhs_evals += solution.num_rhs_evals;
  }
  return RainshaftSolution(states, dvars, num_rhs_evals);
}
