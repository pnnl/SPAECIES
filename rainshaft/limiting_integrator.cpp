#include "limiting_integrator.hpp"
#include <algorithm>

LimitingIntegrator::LimitingIntegrator(const RainshaftConstants &constants,
                                       const RainshaftIntegrator& sub_integrator)
  : constants(constants), integrator(sub_integrator) {
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
    // If qr has to be corrected, reduce q and increase T to allow for better
    // water and energy conservation.
    if (qr[i] < 0.) {
      // Note that for the "original" P3 method applied to the rainshaft model,
      // which this class was designed for, the only way to get significant
      // negative qr is from excessive evaporation (a source of q), so this
      // should never drive q negative. Notably, this is only true so long as
      // the CFL condition is not violated in the sedimentation of qr (which
      // should not be a problem so long as sed_cfl_integrator is used).
      q[i] += qr[i];
      t[i] -= qr[i] * constants.latvap / constants.cp;
      qr[i] = 0.;
    }
    // nr is not conserved by local processes; no special logic required here.
    nr[i] = std::max(nr[i], 0.);
    // Note: This t limiter is not in the original P3 code, which may silently
    // accept negative t or crash on negative t, depending on treatment of
    // floating point exceptions and where the negative t occurs.
    t[i] = std::max(t[i], 0.);
    // In the rainshaft model, large negative q values should not be produced
    // since there is no "sink" except nudging (which should only be used with
    // time scales >> dt), but the below will violate energy/mass conservation
    // if more processes are added without adapting more of P3's limiter logic.
    q[i] = std::max(q[i], 0.);
  }
  return RainshaftSolution({limited_state}, solution.num_rhs_evals);
}
