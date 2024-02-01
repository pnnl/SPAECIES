#include "rainshaft_sum_process.hpp"
#include <iostream>

SumProcess::SumProcess(const std::vector<const RainshaftProcess *>& processes)
  : nsub(processes.size()), sub_processes(processes) {
}

RainshaftTendency SumProcess::calc_tend(const RainshaftConstants& constants,
                                        const RainshaftGrid& grid,
                                        const RainshaftState& state,
                                        const RainshaftDerivedVars& dvars) const {
  std::vector<double> t_tend(grid.nlev, 0.), q_tend(grid.nlev, 0.);
  std::vector<double> nr_tend(grid.nlev, 0.), qr_tend(grid.nlev, 0.);
  // Iterate through sub-processes and add up all their tendencies.
  for (std::size_t is = 0; is != nsub; ++is) {
    RainshaftTendency sub_tend = sub_processes[is]->calc_tend(constants, grid, state, dvars);
    for (std::size_t il = 0; il != grid.nlev; ++il) {
      t_tend[il] += sub_tend.t_tend[il];
      q_tend[il] += sub_tend.q_tend[il];
      nr_tend[il] += sub_tend.nr_tend[il];
      qr_tend[il] += sub_tend.qr_tend[il];
    }
  }
  return RainshaftTendency(t_tend, q_tend, nr_tend, qr_tend);
}


RainshaftTendencyJac SumProcess::calc_tend_jac(const RainshaftConstants& constants,
                                        const RainshaftGrid& grid,
                                        const RainshaftState& state,
                                        const RainshaftDerivedVars& dvars) const {
  double* t_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* q_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* nr_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};
  double* qr_tend_jac = new double[4*grid.nlev * 4*grid.nlev] {0};

  // Iterate through sub-processes and add up all their tendencies.
  for (std::size_t is = 0; is != nsub; ++is) {
    RainshaftTendencyJac sub_tend_jac = sub_processes[is]->calc_tend_jac(constants, grid, state, dvars);
    for (std::size_t il = 0; il != 4*grid.nlev; ++il) {
      for (std::size_t jl = 0; jl != 4*grid.nlev; ++jl) {
        t_tend_jac[il*4*grid.nlev + jl] += sub_tend_jac.t_tend_jac[il*4*grid.nlev + jl];
        q_tend_jac[il*4*grid.nlev + jl] += sub_tend_jac.q_tend_jac[il*4*grid.nlev + jl];
        nr_tend_jac[il*4*grid.nlev + jl] += sub_tend_jac.nr_tend_jac[il*4*grid.nlev + jl];
        qr_tend_jac[il*4*grid.nlev + jl] += sub_tend_jac.qr_tend_jac[il*4*grid.nlev + jl];
      }
    }
  }

  // int nz = 4*grid.nlev;
  // for (int i = nz/2; i != nz; ++i) {
  //   for (int j = nz/2; j != nz; ++j) {
  //     std::cout << qr_tend_jac[i*nz + j] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // std::cout << std::endl;

  return RainshaftTendencyJac(t_tend_jac, q_tend_jac, nr_tend_jac, qr_tend_jac);
}
