#ifndef IMEX_INTEGRATOR_HPP
#define IMEX_INTEGRATOR_HPP

#include <string>
#include <optional>
#include "sundials_integrator.hpp"

class IMEXIntegrator : public SundialsIntegrator<2>
{
private:
  const double dt;
  const int order;
  const double rel_tol;
  const bool postprocess;
  const std::optional<std::string> jacobian_file;

public:
  IMEXIntegrator(const RainshaftConstants &constants,
                     const RainshaftGrid &grid,
                     const RainshaftProcess *const process_imp,
                     const RainshaftProcess *const process_exp,
                     const VarDescList& state_descs,
                     const VarDescList& tend_descs,
                     const double dt = 0,
                     const int order = 4,
                     const double rel_tol = 1.e-4,
                     const bool postprocess = false,
                     const int steps_per_output = -1,
                     const std::optional<std::string> jacobian_file = std::nullopt);

  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const StateConst &initial_state,
                              int& error_flag) const;
};

#endif // IMEX_INTEGRATOR_HPP
