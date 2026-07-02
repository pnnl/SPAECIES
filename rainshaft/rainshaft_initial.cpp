#include "rainshaft_initial.hpp"

bool warm_adiabatic_initial_condition(const RainshaftConstants &constants,
                                      const RainshaftGrid &grid,
                                      const SaturationFormulae &sat_form,
                                      double srf_temp,
                                      double lapse_rate,
                                      double rel_hum_init,
                                      State& initial_state) {
    // Coming up with an initial condition for t and q is slightly tricky, because
  // we only have an implicit relationship between t, q, and dz. But since the
  // effect of q on layer height is not large, start by ignoring it, in which
  // case we do have an explicit relationship between t and dz.
  std::size_t nlev = grid.nlev;
  VarMut t = initial_state.get_variable("T").value();
  VarMut q = initial_state.get_variable("q").value();
  VarMut nr = initial_state.get_variable("nr").value();
  VarMut qr = initial_state.get_variable("qr").value();
  for (std::size_t i = 0; i != nlev; ++i) {
    nr[i] = 0.;
    qr[i] = 0.;
  }
  double rog = constants.rdry / constants.g;
  std::vector<double> z_int(nlev+1, 0.);
  for (int il = nlev - 1; il != -1; --il) {
    double pdel = grid.p_int[il+1]-grid.p_int[il];
    t[il] = (srf_temp - lapse_rate*z_int[il+1])
      / (1. + 0.5 * lapse_rate * rog * pdel/ grid.p_mid[il]);
    z_int[il] = rog * t[il] * pdel / grid.p_mid[il];
  }
  // Now iterate until change in dz between successive iterations is small;
  // ten iterations should do it.
  bool converged = false;
  double t_tolerance = 1.e-3; // Kelvin
  for (int it = 0; it != 10 && !converged; ++it) {
    std::vector<double> t_v(nlev, 0.);
    for (std::size_t il = 0; il != nlev; ++il) {
      q[il] = rel_hum_init * sat_form.q_sat_dry(t[il], grid.p_mid[il]);
      // Virtual temperature factor.
      double t_v_fac = 1. + ((1/constants.epsilon_h2o - 1.) * q[il]);
      t_v[il] = t[il] * t_v_fac;
    }
    auto dz = grid.calc_dz(constants, t_v);
    z_int = dz_to_z_int(dz);
    converged = true;
    for (std::size_t il = 0; il!= nlev; ++il) {
      double new_temp = srf_temp - lapse_rate * (z_int[il] - 0.5 * dz[il]);
      if (abs(new_temp - t[il]) > t_tolerance) {
        converged = false;
      }
      t[il] = new_temp;
    }
  }
  return converged;
}
