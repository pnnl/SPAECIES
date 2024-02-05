#include "rainshaft_sum_process.hpp"

SumProcess::SumProcess(const std::vector<const RainshaftProcess *>& processes)
  : nsub(processes.size()), sub_processes(processes) {
}

void SumProcess::calc_tend(const RainshaftConstants& constants,
                           const RainshaftGrid& grid,
                           const spaecies::VariableArrayView<double>& state,
                           const RainshaftDerivedVars& dvars,
                           spaecies::VariableArrayView<double>& tend) const {
  double* tend_ptr = tend.data();
  spaecies::VariableArray<double> sub_tend(tend.var_descs);
  double* sub_ptr = sub_tend.data();
  // Iterate through sub-processes and add up all their tendencies.
  for (std::size_t is = 0; is != nsub; ++is) {
    for (std::size_t i = 0; i != tend.size; ++i) {
      sub_ptr[i] = 0.;
    }
    sub_processes[is]->calc_tend(constants, grid, state, dvars, sub_tend);
    for (std::size_t i = 0; i != tend.size; ++i) {
      tend_ptr[i] += sub_ptr[i];
    }
  }
}
