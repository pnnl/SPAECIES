#ifndef RAINSHAFT_PROCESS_HPP
#define RAINSHAFT_PROCESS_HPP
#include <vector>

#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_state.hpp"
#include "rainshaft_derived_vars.hpp"
#include "rainshaft_tendency.hpp"

class RainshaftProcess
{

public:
  // Calculate tendency from current state.
  virtual RainshaftTendency calc_tend(const RainshaftConstants &constants,
                                      const RainshaftGrid &grid,
                                      const RainshaftState &state,
                                      const RainshaftDerivedVars &dvars) const = 0;

  virtual void calc_tend_jac_prod(const RainshaftConstants &constants,
                                  const RainshaftGrid &grid,
                                  const RainshaftState &state,
                                  const RainshaftDerivedVars &dvars,
                                  const double *const vec,
                                  double *const prod) const = 0;

  // virtual void tend_jacobian_nnz(const RainshaftGrid& grid) const = 0;

  // virtual void tend_jacobian_sparsity(const RainshaftGrid& grid, const int *colptrs, const int *rowvals) const = 0;
};

#endif // RAINSHAFT_PROCESS_HPP
