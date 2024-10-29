#include "limiting_integrator.hpp"
#include <algorithm>

LimitingIntegrator::LimitingIntegrator(const RainshaftIntegrator& sub_integrator)
  : integrator(sub_integrator) {
}

RainshaftSolution LimitingIntegrator::integrate(double initial_time,
                                                double final_time,
                                                const StateConst& initial_state) const {
  RainshaftSolution solution = integrator.integrate(initial_time, final_time, initial_state);
  State limited_state = solution.states.back().deep_copy();
  VarMut t = limited_state.get_variable("T");
  VarMut q = limited_state.get_variable("q");
  VarMut nr = limited_state.get_variable("nr");
  VarMut qr = limited_state.get_variable("qr");
  for (std::size_t i = 0; i != t.size(); ++i) {
    t[i] = std::max(t[i], 0.);
    q[i] = std::max(q[i], 0.);
    nr[i] = std::max(nr[i], 0.);
    qr[i] = std::max(qr[i], 0.);
  }
  return RainshaftSolution({limited_state}, solution.num_rhs_evals);
}
