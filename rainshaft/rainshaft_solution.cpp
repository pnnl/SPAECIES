#include "rainshaft_solution.hpp"
#include <stdexcept>

RainshaftSolution::RainshaftSolution(const std::vector<spaecies::State<double>>& states,
                                     int num_rhs_evals_in)
  : states(states), num_rhs_evals(num_rhs_evals_in) {
}

RainshaftSolution::RainshaftSolution(const std::vector<spaecies::State<double>>&& states,
                  int num_rhs_evals_in)
  : states(std::move(states)), num_rhs_evals(num_rhs_evals_in) {
}

void RainshaftSolution::pop_back() {
  this->states.pop_back();
}

void RainshaftSolution::move_last_to_other(RainshaftSolution& other) {
  // Avoids copying the VarDescPtr objects, but might be premature optimization...
  other.states.emplace_back(std::vector<spaecies::VarDescPtr>(), nullptr);
  std::swap(other.states.back(), this->states.back());
  this->states.pop_back();
}
