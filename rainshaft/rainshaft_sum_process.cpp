#include "rainshaft_sum_process.hpp"

SumProcess::SumProcess(const std::vector<const RainshaftProcess *> &processes)
    : sub_processes(processes)
{
}

void SumProcess::calc_tend(const RainshaftConstants& constants,
                           const RainshaftGrid& grid,
                           const StateConst& state,
                           const RainshaftDerivedVars& dvars,
                           Tendency& tend) const {

  for (const RainshaftProcess * const p : sub_processes)
  {
    p->calc_tend(constants, grid, state, dvars, tend);
  }
}

void SumProcess::calc_tend_jac(const RainshaftConstants &constants,
                               const RainshaftGrid &grid,
                               const StateConst& state,
                               const RainshaftDerivedVars &dvars,
                               Matrix jac) const
{
  for (const RainshaftProcess * const p : sub_processes)
  {
    p->calc_tend_jac(constants, grid, state, dvars, jac);
  }
}

std::set<std::string> SumProcess::get_required_vars() const
{
  std::set<std::string> vars;
  for (const RainshaftProcess * const p : sub_processes)
  {
    const std::set<std::string> process_vars = p->get_required_vars();
    vars.insert(process_vars.begin(), process_vars.end());
  }
  return vars;
}

std::set<std::string> SumProcess::get_optional_vars() const
{
  std::set<std::string> vars;
  for (const RainshaftProcess * const p : sub_processes)
  {
    const std::set<std::string> process_vars = p->get_optional_vars();
    vars.insert(process_vars.begin(), process_vars.end());
  }
  return vars;
}
