#ifndef RAINSHAFT_SOLUTION_HPP
#define RAINSHAFT_SOLUTION_HPP
#include <vector>

#include "spaecies.hpp"

#include "rainshaft_derived_vars.hpp"

class RainshaftSolution {

public:

  // Constructor from state arrays.
  RainshaftSolution(const std::vector<spaecies::State<const double>>& states,
                    int num_rhs_evals_in);

  // Move constructor from state arrays.
  RainshaftSolution(const std::vector<spaecies::State<const double>>&& states,
                    int num_rhs_evals_in);

  // State vector.
  std::vector<spaecies::State<const double>> states;
  // Number of evaluations of RHS function.
  int num_rhs_evals;

  void pop_back();

  void move_last_to_other(RainshaftSolution& other);
};

#endif // RAINSHAFT_SOLUTION_HPP
