#ifndef SUM_PROCESS_HPP
#define SUM_PROCESS_HPP
#include <cstddef>
#include <array>

#include "rainshaft_process.hpp"

template <std::size_t N>
class SumProcess : public RainshaftProcess
{

public:
  // Constructor from other processes.
  template<typename ... P>
  constexpr SumProcess(P&& ... processes) noexcept
      : sub_processes{std::forward<P>(processes)...} {}

  // Calculate tendency from current state.
  constexpr RainshaftTendency calc_tend(const RainshaftConstants &constants,
                              const RainshaftGrid &grid,
                              const RainshaftState &state,
                              const RainshaftDerivedVars &dvars) const
  {
    std::vector<double> t_tend(grid.nlev, 0.), q_tend(grid.nlev, 0.);
    std::vector<double> nr_tend(grid.nlev, 0.), qr_tend(grid.nlev, 0.);
    // Iterate through sub-processes and add up all their tendencies.
    for (const auto process : sub_processes)
    {
      RainshaftTendency sub_tend = process->calc_tend(constants, grid, state, dvars);
      for (std::size_t il = 0; il != grid.nlev; ++il)
      {
        t_tend[il] += sub_tend.t_tend[il];
        q_tend[il] += sub_tend.q_tend[il];
        nr_tend[il] += sub_tend.nr_tend[il];
        qr_tend[il] += sub_tend.qr_tend[il];
      }
    }
    return RainshaftTendency(t_tend, q_tend, nr_tend, qr_tend);
  }

private:
  // Sub-processes that this process is a sum of.
  const std::array<const RainshaftProcess *, N> sub_processes;
};

// Allow template size deduction
template<typename ... P>
SumProcess(P...) -> SumProcess<sizeof...(P)>;

#endif // SUM_PROCESS_HPP
