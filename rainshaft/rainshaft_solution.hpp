#ifndef RAINSHAFT_SOLUTION_HPP
#define RAINSHAFT_SOLUTION_HPP
#include <vector>

#include "rainshaft_state.hpp"
#include "rainshaft_derived_vars.hpp"

class RainshaftSolution {

public:

  // Constructor from state and derived variable vectors.
  RainshaftSolution(const std::vector<RainshaftState>& state_vec,
                    const std::vector<RainshaftDerivedVars>& dvar_vec);

  // State vector.
  const std::vector<RainshaftState> states;
  // Derived variables vector.
  const std::vector<RainshaftDerivedVars> dvars;

};

#endif // RAINSHAFT_SOLUTION_HPP
