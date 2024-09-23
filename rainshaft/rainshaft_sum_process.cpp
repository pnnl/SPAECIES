#include "rainshaft_sum_process.hpp"

SumProcess::SumProcess(const std::vector<const RainshaftProcess *> &processes)
    : nsub(processes.size()), sub_processes(processes)
{
}

void SumProcess::calc_tend(const RainshaftConstants& constants,
                           const RainshaftGrid& grid,
                           const StateConst& state,
                           const RainshaftDerivedVars& dvars,
                           const Tendency& tend) const {
  // TODO: Why is tend a const reference if it will be modified?
  double* tend_ptr = tend.data();
  Tendency sub_tend(tend.var_descs());
  double* sub_ptr = sub_tend.data();
  // Iterate through sub-processes and add up all their tendencies.
  for (std::size_t is = 0; is != nsub; ++is) {
    // TODO: not needed?
    std::fill_n(sub_ptr, sub_tend.size(), 0.);
    sub_processes[is]->calc_tend(constants, grid, state, dvars, sub_tend);

    // TODO: std::copy_n
    for (std::size_t i = 0; i != tend.size(); ++i) {
      tend_ptr[i] += sub_ptr[i];
    }
  }
}

void SumProcess::calc_tend_jac_prod(const RainshaftConstants &constants,
                                    const RainshaftGrid &grid,
                                    const StateConst& state,
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
                               const StateConst& state,
                               const RainshaftDerivedVars &dvars,
                               Matrix jac) const
{
  for (const auto &p : sub_processes)
  {
    p->calc_tend_jac(constants, grid, state, dvars, jac);
  }
}
