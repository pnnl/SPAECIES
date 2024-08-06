#include "rainshaft_sum_process.hpp"
#include <iostream>

SumProcess::SumProcess(const std::vector<const RainshaftProcess *> &processes)
    : nsub(processes.size()), sub_processes(processes)
{
}

RainshaftTendency SumProcess::calc_tend(const RainshaftConstants &constants,
                                        const RainshaftGrid &grid,
                                        const RainshaftState &state,
                                        const RainshaftDerivedVars &dvars) const
{
  std::vector<double> t_tend(grid.nlev, 0.), q_tend(grid.nlev, 0.);
  std::vector<double> nr_tend(grid.nlev, 0.), qr_tend(grid.nlev, 0.);
  // Iterate through sub-processes and add up all their tendencies.
  for (std::size_t is = 0; is != nsub; ++is)
  {
    RainshaftTendency sub_tend = sub_processes[is]->calc_tend(constants, grid, state, dvars);
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

void SumProcess::calc_tend_jac_prod(const RainshaftConstants &constants,
                                    const RainshaftGrid &grid,
                                    const RainshaftState &state,
                                    const RainshaftDerivedVars &dvars,
                                    const double *const vec,
                                    double *const prod) const
{
  for (const auto &p : sub_processes)
  {
    p->calc_tend_jac_prod(constants, grid, state, dvars, vec, prod);
  }
}

void SumProcess::calc_tend_jac(const RainshaftConstants &constants,
                               const RainshaftGrid &grid,
                               const RainshaftState &state,
                               const RainshaftDerivedVars &dvars,
                               Matrix jac) const
{
  for (const auto &p : sub_processes)
  {
    p->calc_tend_jac(constants, grid, state, dvars, jac);
  }
}
