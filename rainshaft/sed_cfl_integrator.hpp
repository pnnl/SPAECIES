#ifndef SED_CFL_INTEGRATOR_HPP
#define SED_CFL_INTEGRATOR_HPP

#include "spaecies.hpp"

#include "rainshaft_integrator.hpp"
#include "sedimentation.hpp"

class SedCflIntegrator : public RainshaftIntegrator {

public:

  SedCflIntegrator(const RainshaftConstants* constants,
                   const RainshaftGrid* grid,
                   const Sedimentation *sedimentation,
                   const RainshaftIntegrator *integrator);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const spaecies::VariableArray<double>& initial_state) const;

private:

  const RainshaftConstants* constants;

  const RainshaftGrid* grid;

  const Sedimentation *sed;

  const RainshaftIntegrator *sed_int;

};

#endif // SED_CFL_INTEGRATOR_HPP
