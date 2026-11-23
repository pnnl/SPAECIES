#ifndef RAINSHAFT_SOLUTION_HPP
#define RAINSHAFT_SOLUTION_HPP
#include <vector>

#include "rainshaft_derived_vars.hpp"
#include "rainshaft_types.hpp"

class RainshaftSolution {

public:

  // Constructor from state arrays.
  RainshaftSolution(const std::vector<StateConst>& states,
                    const std::vector<double>& times,
                    int num_rhs_evals_in);

  // Move constructor from state arrays.
  RainshaftSolution(std::vector<StateConst>&& states,
                    std::vector<double>&& times,
                    int num_rhs_evals_in);

  // State vector.
  std::vector<StateConst> states;
  // Output time corresponding to each state.
  std::vector<double> times;
  // Number of evaluations of RHS function.
  int num_rhs_evals;

  void pop_back();

  void move_last_to_other(RainshaftSolution& other);
};

#endif // RAINSHAFT_SOLUTION_HPP
