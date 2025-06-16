#ifndef SED_CFL_INTEGRATOR_HPP
#define SED_CFL_INTEGRATOR_HPP

#include "rainshaft_integrator.hpp"
#include "sedimentation.hpp"

class SedCflIntegrator : public RainshaftIntegrator {

public:

  SedCflIntegrator(const RainshaftConstants& constants,
                   const RainshaftGrid& grid,
                   const VarDescList& tend_descs,
                   const Sedimentation& sedimentation);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const StateConst& initial_state) const;

private:

  const RainshaftConstants& constants;

  const RainshaftGrid& grid;

  const VarDescList tend_descs;

  const Sedimentation &sed;

};

#endif // SED_CFL_INTEGRATOR_HPP
