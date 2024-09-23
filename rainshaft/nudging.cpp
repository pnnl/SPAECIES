#include "nudging.hpp"
#include <cstddef>
#include <stdexcept>

Nudging::Nudging(double time_scale_in, const std::vector<double> &t, const std::vector<double> &q)
    : time_scale(time_scale_in), t0(t), q0(q)
{
  // Check vector sizes.
  // SPS: Should actually print out mismatched dimensions on failure.
  if (t.size() != q.size())
  {
    throw std::invalid_argument("t and q differ in size");
  }
}

void Nudging::calc_tend(const RainshaftConstants& constants,
                        const RainshaftGrid& grid,
                        const StateConst& state,
                        const RainshaftDerivedVars& dvars,
                        const Tendency& tend) const {
  VarConst t = state.get_variable("T");
  VarConst q = state.get_variable("q");
  VarMut t_tend = tend.get_variable("T_tend");
  VarMut q_tend = tend.get_variable("q_tend");
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    t_tend[il] = (t0[il] - t[il]) / time_scale;
    q_tend[il] = (q0[il] - q[il]) / time_scale;
  }
}

void Nudging::calc_tend_jac_prod(const RainshaftConstants &constants,
                                 const RainshaftGrid &grid,
                                 const StateConst& state,
                                 const RainshaftDerivedVars &dvars,
                                 const double *const vec,
                                 double *const prod) const
{
  for (std::size_t il = 0; il != 2 * grid.nlev; ++il)
  {
    prod[il] -= vec[il] / time_scale;
  }
}

void Nudging::calc_tend_jac(const RainshaftConstants &constants,
                            const RainshaftGrid &grid,
                            const StateConst& state,
                            const RainshaftDerivedVars &dvars,
                            Matrix jac) const
{
  for (std::size_t il = 0; il != 2 * grid.nlev; ++il)
  {
    jac(il, il) -= 1.0 / time_scale;
  }
}
