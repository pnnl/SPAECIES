#include "self_collision.hpp"
#include <cstddef>
#include <cmath>
using std::min, std::cbrt, std::exp;

RainshaftTendency SelfCollision::calc_tend(const RainshaftConstants& constants,
                                           const RainshaftGrid& grid,
                                           const RainshaftState& state,
                                           const RainshaftDerivedVars& dvars) const {
  std::vector<double> t_tend(grid.nlev, 0.), q_tend(grid.nlev, 0.);
  std::vector<double> nr_tend(grid.nlev, 0.), qr_tend(grid.nlev, 0.);
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    double b = breakup_fac(constants, state.nr[il], state.qr[il]);
    nr_tend[il] = -5.78 * b * state.nr[il] * state.qr[il] * dvars.rho_dry[il];
  }
  return RainshaftTendency(t_tend, q_tend, nr_tend, qr_tend);
}

// placeholder for jacobian of self collision!
RainshaftTendencyJac SelfCollision::calc_tend_jac(const RainshaftConstants& constants,
                                           const RainshaftGrid& grid,
                                           const RainshaftState& state,
                                           const RainshaftDerivedVars& dvars) const {
  double* t_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* q_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* nr_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* qr_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};

  // d(rsc)/dT
  for (std::size_t il = 2*grid.nlev; il < 3*grid.nlev; ++il) {
    for (std::size_t jl = 0; jl < grid.nlev; ++jl) {
      if (il == jl+2*grid.nlev) {
        double b = breakup_fac(constants, state.nr[jl], state.qr[jl]);
        t_tend_jac[jl*4*grid.nlev + il] = -5.78 * b * state.nr[jl] * state.qr[jl] * 
                                            grid.p_mid[jl]/(constants.rdry * (1.0 + state.q[jl]/constants.epsilon_h2o)) * -1.0/pow(state.t[jl], 2);
      }
    }
  }

  // d(rsc)/dq
  for (std::size_t il = 2*grid.nlev; il < 3*grid.nlev; ++il) {
    for (std::size_t jl = grid.nlev; jl < 2*grid.nlev; ++jl) {
      if (il == jl+grid.nlev) {
        double b = breakup_fac(constants, state.nr[jl-grid.nlev], state.qr[jl-grid.nlev]);
        q_tend_jac[jl*4*grid.nlev + il] = -5.78 * b * state.nr[jl-grid.nlev] * state.qr[jl-grid.nlev] * 
                                            grid.p_mid[jl-grid.nlev]/(constants.rdry*state.t[jl-grid.nlev]*constants.epsilon_h2o) * -1.0/pow(1.0 + state.q[jl-grid.nlev]/constants.epsilon_h2o, 2);
      }
    }
  }

  // d(rsc)/d(nr)
  for (std::size_t il = 2*grid.nlev; il < 3*grid.nlev; ++il) {
    for (std::size_t jl = 2*grid.nlev; jl < 3*grid.nlev; ++jl) {
      if (il == jl) {
        double b = breakup_fac(constants, state.nr[jl-2*grid.nlev], state.qr[jl-2*grid.nlev]);
        double db;

        double mean_mass_diam = std::cbrt(state.qr[jl-2*grid.nlev] / (constants.pi * constants.rhow * state.nr[jl-2*grid.nlev]));
        double F = 2. - std::exp(2300. * (2.8e-4 - mean_mass_diam));
        if (1.0 <= F) {
          db = 0.0;
        } else {
          if (dvars.lambdar[jl-2*grid.nlev] == 0.0) {
            db = - std::exp(2300. * (2.8e-4 - 1.0/(1.e-5))) * -2300. / (3.0*(1.e-5));
          } else {
            db = - std::exp(2300. * (2.8e-4 - 1.0/dvars.lambdar[jl-2*grid.nlev])) * -2300. / (3.0*dvars.lambdar[jl-2*grid.nlev]*state.nr[jl-2*grid.nlev]);
          }
        }
        nr_tend_jac[jl*4*grid.nlev + il] = -5.78 * state.qr[jl-2*grid.nlev] * dvars.rho_dry[jl-2*grid.nlev] * (db * state.nr[jl-2*grid.nlev] + b);
                                          
      }
    }
  }


  // d(rsc)/d(qr)
  for (std::size_t il = 2*grid.nlev; il < 3*grid.nlev; ++il) {
    for (std::size_t jl = 3*grid.nlev; jl < 4*grid.nlev; ++jl) {
      if (il == jl-grid.nlev) {
        double b = breakup_fac(constants, state.nr[jl-3*grid.nlev], state.qr[jl-3*grid.nlev]);
        double db;

        double mean_mass_diam = std::cbrt(state.qr[jl-3*grid.nlev] / (constants.pi * constants.rhow * state.nr[jl-3*grid.nlev]));
        double F = 2. - std::exp(2300. * (2.8e-4 - mean_mass_diam));
        if (1.0 <= F) {
          db = 0.0;
        } else {
          if (dvars.lambdar[jl-3*grid.nlev] == 0.0) {
            db = - std::exp(2300. * (2.8e-4 - 1.0/(1.e-5))) * -2300. / (-3.0*(1.e-5));
          } else {
            db = - std::exp(2300. * (2.8e-4 - 1.0/dvars.lambdar[jl-3*grid.nlev])) * -2300. / (-3.0*dvars.lambdar[jl-3*grid.nlev]*state.qr[jl-3*grid.nlev]);
          }
        }
        qr_tend_jac[jl*4*grid.nlev + il] = -5.78 * state.nr[jl-3*grid.nlev] * dvars.rho_dry[jl-3*grid.nlev] * (db * state.qr[jl-3*grid.nlev] + b);
                                          
      }
    }
  }


  return RainshaftTendencyJac(t_tend_jac, q_tend_jac, nr_tend_jac, qr_tend_jac);
}


double SelfCollision::breakup_fac(const RainshaftConstants& constants,
                                               double nr, double qr) const {
  // SPS: Should put in an assert that nr is sufficiently large?
  // Diameter of a particle of mean mass (m).
  double mean_mass_diam = std::cbrt(qr / (constants.pi * constants.rhow * nr));
  // Returned factor (dimensionless).
  double b = min(1., 2. - std::exp(2300. * (2.8e-4 - mean_mass_diam)));
  // double b = 2. - std::exp(2300. * (2.8e-4 - mean_mass_diam));
  return b;
}
