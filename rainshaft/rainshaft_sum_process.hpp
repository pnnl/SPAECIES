#ifndef SUM_PROCESS_HPP
#define SUM_PROCESS_HPP
#include <cstddef>
#include <vector>

#include "rainshaft_process.hpp"

class SumProcess : public RainshaftProcess {

public:

  // Constructor from other processes.
  // SPS: Need to consider other pointer types for proper memory safety.
  SumProcess(const std::vector<const RainshaftProcess *>& processes);

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

private:

  // Sub-processes that this process is a sum of.
  const std::vector<const RainshaftProcess *> sub_processes;

};

#endif // SUM_PROCESS_HPP
