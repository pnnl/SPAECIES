#ifndef CLOUD_SEDIMENTATION_HPP
#define CLOUD_SEDIMENTATION_HPP
#include <vector>
#include <optional>
#include <tuple>
#include <boost/math/special_functions/gamma.hpp>
using boost::math::tgamma;

#include "rainshaft_process.hpp"

class CloudSedimentation : public RainshaftProcess 
{
private:
  double calc_moment_flux(const double moment,
                          const double rho_dry,
                          const double v) const;

  double calc_moment_tend(const double dz,
                          const double rho_dry,
                          const double flux_bot,
                          const double flux_top) const;

  double calc_muc(double nc, double rho_dry) const;

  double calc_lambdac(const RainshaftConstants& constants, double nc, double qc, double muc) const;

public:

  CloudSedimentation();

  // Calculate tendency from current state.
  void calc_tend(const RainshaftConstants& constants,
                 const RainshaftGrid& grid,
                 const StateConst& state,
                 const RainshaftDerivedVars& dvars,
                 Tendency& tend) const;

  void calc_tend_jac(const RainshaftConstants &constants,
                     const RainshaftGrid &grid,
                     const StateConst& state,
                     const RainshaftDerivedVars &dvars,
                     Matrix jac) const;

  // For a given size distribution+temperature, what are the cloud number and mass fall speeds?
  std::tuple<double, double> cloud_fall_speeds(const RainshaftConstants &constants, 
                                               double t, double lambdac, double muc) const;

};

#endif // CLOUD_SEDIMENTATION_HPP
