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

RainshaftTendency Nudging::calc_tend(const RainshaftConstants& constants,
                                     const RainshaftGrid& grid,
                                     const RainshaftState& state,
                                     const RainshaftDerivedVars& dvars) const {
  std::vector<double> t_tend(grid.nlev, 0.), q_tend(grid.nlev, 0.);
  std::vector<double> nr_tend(grid.nlev, 0.), qr_tend(grid.nlev, 0.);
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    t_tend[il] = (t0[il] - state.t[il]) / time_scale;
    q_tend[il] = (q0[il] - state.q[il]) / time_scale;
  }
  return RainshaftTendency(t_tend, q_tend, nr_tend, qr_tend);
}

// placeholder for Jacobian!
RainshaftTendencyJac Nudging::calc_tend_jac(const RainshaftConstants& constants,
                                     const RainshaftGrid& grid,
                                     const RainshaftState& state,
                                     const RainshaftDerivedVars& dvars) const {
  double* t_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* q_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* nr_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* qr_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};

  // d(nudge_T)/dT
  for (std::size_t il = 0; il < grid.nlev; ++il) {
    for (std::size_t jl = 0; jl < grid.nlev; ++jl) {
      if (il == jl) {
        t_tend_jac[jl*4*grid.nlev + il] = -1.0 / time_scale;
      }
    }
  }

  // d(nudge_q)/dT
  for (std::size_t il = grid.nlev; il < 2*grid.nlev; ++il) {
    for (std::size_t jl = grid.nlev; jl < 2*grid.nlev; ++jl) {
      if (il == jl) {
        q_tend_jac[jl*4*grid.nlev + il] = -1.0 / time_scale;
      }
    }
  }

  return RainshaftTendencyJac(t_tend_jac, q_tend_jac, nr_tend_jac, qr_tend_jac);
}
