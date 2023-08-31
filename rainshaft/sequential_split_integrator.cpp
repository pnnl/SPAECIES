#include "sequential_split_integrator.hpp"

SequentialSplitIntegrator::SequentialSplitIntegrator(const std::vector<const RainshaftIntegrator *>& sub_integrators)
  : integrators(sub_integrators) {
}

RainshaftSolution SequentialSplitIntegrator::integrate(double initial_time,
                                                       double final_time,
                                                       const RainshaftState& initial_state) const {
  std::vector<RainshaftState> states{initial_state};
  std::vector<RainshaftDerivedVars> dvars;
  // This is a very inelegant way to keep the latest solution in scope,
  // since solutions have const members and can't be copied over.
  std::vector<RainshaftSolution> solutions;
  int num_rhs_evals = 0;
  RainshaftSolution first_solution = integrators[0]->integrate(initial_time, final_time, initial_state);
  // Handle initial dvars, since this class can't directly calculate it.
  // SPS: Need a different strategy if there are no sub-integrators?
  dvars.push_back(first_solution.dvars[0]);
  num_rhs_evals += first_solution.num_rhs_evals;
  solutions.push_back(first_solution);
  // Loop over all integrators other than the first.
  for (std::size_t ii = 1; ii != integrators.size(); ++ii) {
    RainshaftSolution solution = integrators[ii]->integrate(initial_time, final_time,
                                                            solutions.back().states.back());
    num_rhs_evals += solution.num_rhs_evals;
    solutions.pop_back();
    solutions.push_back(solution);
  }
  states.push_back(solutions.back().states.back());
  dvars.push_back(solutions.back().dvars.back());
  return RainshaftSolution(states, dvars, num_rhs_evals);
}
