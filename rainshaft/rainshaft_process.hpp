#ifndef RAINSHAFT_PROCESS_HPP
#define RAINSHAFT_PROCESS_HPP
#include <vector>
#include <functional>

#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_derived_vars.hpp"
#include "rainshaft_types.hpp"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_matrix.h"

class RainshaftProcess
{

public:
  using Matrix = std::function<double &(std::size_t, std::size_t)>;

  // Calculate tendency from current state.
  virtual void calc_tend(const RainshaftConstants &constants,
                                      const RainshaftGrid &grid,
                                      const StateConst& state,
                                      const RainshaftDerivedVars &dvars,
                                      const Tendency& tend) const = 0;

  virtual void calc_tend_jac_prod(const RainshaftConstants &constants,
                                  const RainshaftGrid &grid,
                                  const StateConst& state,
                                  const RainshaftDerivedVars &dvars,
                                  const double *const vec,
                                  double *const prod) const = 0;

  virtual void calc_tend_jac(const RainshaftConstants &constants,
                             const RainshaftGrid &grid,
                             const StateConst& state,
                             const RainshaftDerivedVars &dvars,
                             Matrix jac) const = 0;
};

#endif // RAINSHAFT_PROCESS_HPP
