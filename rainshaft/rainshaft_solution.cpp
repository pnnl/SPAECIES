#include "rainshaft_solution.hpp"
#include <stdexcept>

RainshaftSolution::RainshaftSolution(const std::vector<StateConst>& states,
                                     const std::vector<double>& times,
                                     int num_rhs_evals_in)
  : states(states), times(times), num_rhs_evals(num_rhs_evals_in) {
}

RainshaftSolution::RainshaftSolution(std::vector<StateConst>&& states,
                                     std::vector<double>&& times,
                                     int num_rhs_evals_in)
  : states(std::move(states)), times(std::move(times)),
    num_rhs_evals(num_rhs_evals_in) {
}

void RainshaftSolution::pop_back() {
  this->states.pop_back();
  this->times.pop_back();
}

void RainshaftSolution::move_last_to_other(RainshaftSolution& other) {
  other.states.emplace_back(std::move(this->states.back()));
  other.times.emplace_back(this->times.back());
  this->states.pop_back();
  this->times.pop_back();
}
