#include "rainshaft_solution.hpp"
#include <stdexcept>

RainshaftSolution::RainshaftSolution(const std::vector<spaecies::VariableArray<double>>& state_vec,
                                     const std::vector<RainshaftDerivedVars>& dvar_vec,
                                     int num_rhs_evals_in)
  : states(state_vec), dvars(dvar_vec), num_rhs_evals(num_rhs_evals_in) {
  // Check vector sizes.
  // SPS: Should actually print out mismatched dimensions on failure.
  if (state_vec.size() != dvar_vec.size()) {
    throw std::invalid_argument("state_vec and dvar_vec differ in size");
  }
}
