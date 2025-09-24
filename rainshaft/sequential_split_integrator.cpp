#include "sequential_split_integrator.hpp"

SequentialSplitIntegrator::SequentialSplitIntegrator(const std::vector<const RainshaftIntegrator *>& sub_integrators)
  : integrators(sub_integrators) {
}

RainshaftSolution SequentialSplitIntegrator::integrate(double initial_time,
                                                       double final_time,
                                                       const StateConst& initial_state,
                                                       int& error_flag) const {
  int num_rhs_evals = 0;
  RainshaftSolution solution = integrators[0]->integrate(initial_time, final_time, initial_state, error_flag);
  num_rhs_evals += solution.num_rhs_evals;
  // Loop over all integrators other than the first.
  for (std::size_t ii = 1; ii != integrators.size(); ++ii) {
    RainshaftSolution split_solution = integrators[ii]->integrate(initial_time, final_time,
                                                                  solution.states.back(), error_flag);
    num_rhs_evals += split_solution.num_rhs_evals;
    solution.pop_back();
    split_solution.move_last_to_other(solution);
  }
  solution.num_rhs_evals = num_rhs_evals;
  return solution;
}
