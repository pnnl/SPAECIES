#ifndef NUDGING_HPP
#define NUDGING_HPP
#include <vector>

#include "rainshaft_process.hpp"

class Nudging : public RainshaftProcess {

public:

  Nudging(double time_scale_in, const std::vector<double>& t,
          const std::vector<double>& q);

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

  std::set<std::string> get_required_vars() const;
  std::set<std::string> get_optional_vars() const;

  double time_scale;
  const std::vector<double> t0;
  const std::vector<double> q0;

};

#endif // EVAPORATION_HPP
