#ifndef SED_CFL_INTEGRATOR_HPP
#define SED_CFL_INTEGRATOR_HPP

#include "rainshaft_integrator.hpp"
#include "rain_sedimentation.hpp"
#include "size_limiters.hpp"

class SedCflIntegrator : public RainshaftIntegrator {

public:

  SedCflIntegrator(const RainshaftConstants& constants,
                   const RainshaftGrid& grid,
                   const SizeLimiters& size_limiters,
                   const VarDescList& tend_descs,
                   const RainSedimentation& sedimentation,
                   const bool regularize_lambdar);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const StateConst& initial_state,
                              int& error_flag) const;

private:

  const RainshaftConstants& constants;

  const RainshaftGrid& grid;

  const SizeLimiters size_limiters;

  const VarDescList tend_descs;

  const RainSedimentation &sed;

  const bool regularize_lambdar;

};

#endif // SED_CFL_INTEGRATOR_HPP
