#ifndef RAINSHAFT_SOLUTION_HPP
#define RAINSHAFT_SOLUTION_HPP
#include <vector>

#include "spaecies.hpp"

#include "rainshaft_derived_vars.hpp"

class RainshaftSolution {

public:

  // Constructor from state and derived variable vectors.
  RainshaftSolution(const std::vector<spaecies::VariableArray<double>>& state_vec,
                    const std::vector<RainshaftDerivedVars>& dvar_vec,
                    int num_rhs_evals_in);

  // State vector.
  const std::vector<spaecies::VariableArray<double>> states;
  // Derived variables vector.
  const std::vector<RainshaftDerivedVars> dvars;
  // Number of evaluations of RHS function.
  int num_rhs_evals;

};

#endif // RAINSHAFT_SOLUTION_HPP
