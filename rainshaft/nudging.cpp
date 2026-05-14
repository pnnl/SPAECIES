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

void Nudging::calc_tend(const RainshaftConstants&,
                        const RainshaftGrid& grid,
                        const StateConst& state,
                        const RainshaftDerivedVars&,
                        Tendency& tend) const {
  VarConst t = state.get_variable("T");
  VarConst q = state.get_variable("q");
  VarMut t_tend = tend.get_variable("T_tend");
  VarMut q_tend = tend.get_variable("q_tend");
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    t_tend[il] += (t0[il] - t[il]) / time_scale;
    q_tend[il] += (q0[il] - q[il]) / time_scale;
  }
}

void Nudging::calc_tend_jac(const RainshaftConstants &,
                            const RainshaftGrid &grid,
                            const StateConst& state,
                            const RainshaftDerivedVars &,
                            Matrix jac) const
{
  const std::size_t offset_t = state.get_idx("T");
  for (std::size_t il = offset_t; il != offset_t + grid.nlev; ++il)
  {
    jac(il, il) -= 1.0 / time_scale;
  }

  const std::size_t offset_q = state.get_idx("q");
  for (std::size_t il = offset_q; il != offset_q + grid.nlev; ++il)
  {
    jac(il, il) -= 1.0 / time_scale;
  }
}

std::set<std::string> Nudging::get_required_vars() const
{
  return {"T", "q"};
}

std::set<std::string> Nudging::get_optional_vars() const
{
  return {};
}
