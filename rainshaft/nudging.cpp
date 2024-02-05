#include "nudging.hpp"
#include <cstddef>
#include <stdexcept>

Nudging::Nudging(double time_scale_in, const std::vector<double>& t, const std::vector<double>& q)
  : time_scale(time_scale_in), t0(t), q0(q) {
  // Check vector sizes.
  // SPS: Should actually print out mismatched dimensions on failure.
  if (t.size() != q.size()) {
    throw std::invalid_argument("t and q differ in size");
  }
}

void Nudging::calc_tend(const RainshaftConstants& constants,
                        const RainshaftGrid& grid,
                        const spaecies::VariableArrayView<double>& state,
                        const RainshaftDerivedVars& dvars,
                        spaecies::VariableArrayView<double>& tend) const {
  auto t = state.get_variable("T");
  auto q = state.get_variable("q");
  auto t_tend = tend.get_variable("T_tend");
  auto q_tend = tend.get_variable("q_tend");
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    t_tend[il] = (t0[il] - t[il]) / time_scale;
    q_tend[il] = (q0[il] - q[il]) / time_scale;
  }
}
