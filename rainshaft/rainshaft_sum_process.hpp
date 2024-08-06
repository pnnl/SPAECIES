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
  RainshaftTendency calc_tend(const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const RainshaftState& state,
                              const RainshaftDerivedVars& dvars) const;

  void calc_tend_jac_prod(const RainshaftConstants &constants,
                          const RainshaftGrid &grid,
                          const RainshaftState &state,
                          const RainshaftDerivedVars &dvars,
                          const double *const vec,
                          double *const prod) const;

  void calc_tend_jac(const RainshaftConstants &constants,
                             const RainshaftGrid &grid,
                             const RainshaftState &state,
                             const RainshaftDerivedVars &dvars,
                             Matrix jac) const;

private:

  // Number of subprocesses in this process.
  std::size_t nsub;

  // Sub-processes that this process is a sum of.
  std::vector<const RainshaftProcess *> sub_processes;

};

#endif // SUM_PROCESS_HPP
