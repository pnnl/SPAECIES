#include "rainshaft_integrator.hpp"
#include <cstddef>
#include <vector>

RainshaftIntegrator::RainshaftIntegrator(double dt_in) : dt(dt_in) {
}

RainshaftState RainshaftIntegrator::apply_tend(const RainshaftState& state,
                                               const RainshaftTendency& tend) {
  std::vector<double> t, q, nr, qr;
  for (std::size_t i = 0; i != state.t.size(); ++i) {
    t.push_back(state.t[i] + dt*tend.t_tend[i]);
    q.push_back(state.q[i] + dt*tend.q_tend[i]);
    nr.push_back(state.nr[i] + dt*tend.nr_tend[i]);
    qr.push_back(state.qr[i] + dt*tend.qr_tend[i]);
  }
  return RainshaftState(t, q, nr, qr);
}
