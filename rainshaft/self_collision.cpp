#include "self_collision.hpp"
#include "derivatives.hpp"
#include <stdexcept>

SelfCollision::SelfCollision(const bool regularize_lambdar)
  : regularize_lambdar(regularize_lambdar)
{
}

void SelfCollision::calc_tend(const RainshaftConstants& constants,
                              const RainshaftGrid& grid,
                              const StateConst& state,
                              const RainshaftDerivedVars& dvars,
                              Tendency& tend) const {
  VarConst nr = state.get_variable("nr");
  VarConst qr = state.get_variable("qr");
  VarMut nr_tend = tend.get_variable("nr_tend");
  for (std::size_t il = 0; il != grid.nlev; ++il) {
    const auto breakup = breakup_fac(constants, nr[il], qr[il]);
    nr_tend[il] += calc_nr_tend(nr[il], qr[il], breakup, dvars.rho_dry[il]);
  }
}

void SelfCollision::calc_tend_jac(const RainshaftConstants &constants,
                                  const RainshaftGrid &grid,
                                  const StateConst& state,
                                  const RainshaftDerivedVars &dvars,
                                  Matrix jac) const
{
  VarConst T = state.get_variable("T");
  VarConst q = state.get_variable("q");
  VarConst nr = state.get_variable("nr");
  VarConst qr = state.get_variable("qr");

  const std::size_t offset_t = state.get_idx("T");
  const std::size_t offset_q = state.get_idx("q");
  const std::size_t offset_nr = state.get_idx("nr");
  const std::size_t offset_qr = state.get_idx("qr");

  for (std::size_t il = 0; il != grid.nlev; ++il)
  {
    const std::size_t i_t = offset_t + il;
    const std::size_t i_q = offset_q + il;
    const std::size_t i_nr = offset_nr + il;
    const std::size_t i_qr = offset_qr + il;

    const auto breakup = breakup_fac<true>(constants, nr[il], qr[il]);
    const auto rho_dry = dvars.get_rho_dry<true>(constants, T[il], q[il], il);
    const auto nr_tend = calc_nr_tend<true>(nr[il], qr[il], breakup, rho_dry);
    const auto [nr_tend_T, nr_tend_q, nr_tend_nr, nr_tend_qr] = get_grad(nr_tend);

    jac(i_nr, i_t) += nr_tend_T;
    jac(i_nr, i_q) += nr_tend_q;
    jac(i_nr, i_nr) += nr_tend_nr;
    jac(i_nr, i_qr) += nr_tend_qr;
  }
}

std::set<std::string> SelfCollision::get_required_vars() const
{
  return {"T", "q", "nr", "qr"};
}

std::set<std::string> SelfCollision::get_optional_vars() const
{
  return {};
}
