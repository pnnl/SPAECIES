#include "sedimentation.hpp"
#include <cstddef>

RainshaftTendency Sedimentation::calc_tend(const RainshaftConstants& constants,
                                           const RainshaftGrid& grid,
                                           const RainshaftState& state,
                                           const RainshaftDerivedVars& dvars) {
  std::vector<double> t_tend, q_tend, nr_tend, qr_tend;
  t_tend.push_back(0.);
  q_tend.push_back(0.);
  nr_tend.push_back(1.e6 * 2. - 2.*state.nr[0]*dvars.rho[0]);
  qr_tend.push_back(0.003 * 2. - 2.*state.qr[0]*dvars.rho[0]);
  for (std::size_t i = 1; i != grid.nlev; ++i) {
    t_tend.push_back(0.);
    q_tend.push_back(0.);
    nr_tend.push_back(2.*state.nr[i-1]*dvars.rho[i-1]- 2.*state.nr[i]*dvars.rho[i]);
    qr_tend.push_back(2.*state.qr[i-1]*dvars.rho[i-1]- 2.*state.qr[i]*dvars.rho[i]);
  }
  return RainshaftTendency(t_tend, q_tend, nr_tend, qr_tend);
}
