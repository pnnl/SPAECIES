#ifndef SUNDIALS_INTEGRATOR_HPP
#define SUNDIALS_INTEGRATOR_HPP
#include <array>
#include <stdexcept>
#include <limits>
#include "rainshaft_integrator.hpp"
#include "sundials/sundials_context.hpp"
#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "rainshaft_integrator.hpp"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_process.hpp"
#include "rainshaft_solution.hpp"
#include "rainshaft_types.hpp"
#include "variable_array_view.hpp"
#include "size_limiters.hpp"

static_assert(std::is_same_v<sunrealtype, double>, "sunrealtype must be double");

template <int PARTITIONS>
class SundialsIntegrator : public RainshaftIntegrator
{
private:
  using PartitionArray = std::array<const RainshaftProcess *const, PARTITIONS>;

  struct RainshaftUserData
  {
    const RainshaftConstants &constants;
    const RainshaftGrid &grid;
    const SizeLimiters &size_limiters;
    const PartitionArray processes;
    const VarDescList state_descs;
    const VarDescList tend_descs;
    const bool regularize_lambdar;
  };

  template <int PARTITION>
  static int rainshaft_f(sunrealtype, N_Vector y, N_Vector ydot, void *user_data)
  {
    const RainshaftUserData& cast_data = *static_cast<RainshaftUserData *>(user_data);
    const StateConst state = n_vector_to_state(y, cast_data.state_descs);
    const RainshaftDerivedVars dvars = RainshaftDerivedVars(cast_data.constants,
                                                      cast_data.grid,
                                                      state,
                                                      cast_data.regularize_lambdar);
    // Zero out ydot so that we don't have to remember to zero every value in calc_tend.
    N_VConst(0., ydot);
    Tendency tend = n_vector_to_tendency(ydot, cast_data.tend_descs);

    try {
      cast_data.processes[PARTITION]->calc_tend(cast_data.constants,
                                               cast_data.grid,
                                               state, dvars, tend);
    } catch (const std::domain_error &) {
      // Attempt recovery on domain error, e.g., gamma function with negative
      // input
      return 1;
    }
    return 0;
  }

  template <int PARTITION>
  static int rainshaft_jac(sunrealtype, N_Vector y, N_Vector, SUNMatrix Jac, void *user_data, N_Vector, N_Vector, N_Vector)
  {
    SUNMatZero(Jac);
    const RainshaftUserData& cast_data = *static_cast<RainshaftUserData *>(user_data);
    const StateConst state = n_vector_to_state(y, cast_data.state_descs);
    const RainshaftDerivedVars dvars = RainshaftDerivedVars(cast_data.constants,
                                      cast_data.grid,
                                      state,
                                      cast_data.regularize_lambdar);

    RainshaftProcess::Matrix mat = [Jac](const std::size_t i, const std::size_t j) -> sunrealtype & {
      return SM_ELEMENT_D(Jac, i, j);
    };

    try {
      cast_data.processes[PARTITION]->calc_tend_jac(cast_data.constants,
                                                    cast_data.grid,
                                                    state, dvars,
                                                    mat);
    } catch (const std::domain_error &) {
      // Attempt recovery on domain error, e.g., gamma function with negative
      // input
      return 1;
    }
    return 0;
  }

  template <int PARTITION, typename P>
  static int rainshaft_stability(N_Vector y, sunrealtype, sunrealtype *hstab, void *user_data) {
    
    const RainshaftUserData& cast_data = *static_cast<RainshaftUserData *>(user_data);
    const StateConst state = n_vector_to_state(y, cast_data.state_descs);
    const RainshaftDerivedVars dvars = RainshaftDerivedVars(cast_data.constants,
                                      cast_data.grid,
                                      state,
                                      cast_data.regularize_lambdar);
    const auto &process = *static_cast<const P*>(cast_data.processes[PARTITION]);
    *hstab = process.calc_max_step(cast_data.constants, cast_data.grid, dvars);
    return 0;
  }

public:
  SundialsIntegrator(const RainshaftConstants &constants,
                     const RainshaftGrid &grid,
                     const SizeLimiters &size_limiters,
                     PartitionArray processes,
                     const VarDescList& state_descs,
                     const VarDescList& tend_descs,
                     const int steps_per_output,
                     const bool regularize_lambdar)
      : user_data{constants, grid, size_limiters, std::move(processes), state_descs, tend_descs, regularize_lambdar}, 
      steps_per_output(steps_per_output)
  {}

protected:
  const RainshaftUserData user_data;
  const int steps_per_output;
  const sundials::Context sun_ctxt;

  template <int PARTITION>
  auto create_f() const
  {
    return (user_data.processes[PARTITION] == nullptr) ? nullptr : rainshaft_f<PARTITION>;
  }

  template <int PARTITION>
  auto create_jac() const
  {
    return (user_data.processes[PARTITION] == nullptr) ? nullptr : rainshaft_jac<PARTITION>;
  }

  template <int PARTITION, typename P>
  auto create_stability() const
  {
    return (user_data.processes[PARTITION] == nullptr) ? nullptr : rainshaft_stability<PARTITION, P>;
  }

  template <typename E, typename C, typename O, typename Mode>
  RainshaftSolution evolve(E evolveFun,
                           C countFun,
                           O outputFun,
                           void *mem,
                           const double final_time,
                           const N_Vector y,
                           const Mode normal,
                           const Mode one_step,
                           int& error_flag) const
  {
    std::vector<StateConst> states;
    sunrealtype tret = -std::numeric_limits<sunrealtype>::infinity();
    if (steps_per_output > 0)
    {
      for (int i = 0; tret < final_time; i++)
      {
        if (i % steps_per_output == 0) {
          states.emplace_back(n_vector_to_state(y, user_data.state_descs).deep_copy());
          outputFun();
        }
        error_flag = evolveFun(mem, final_time, y, &tret, one_step);
      }
    }
    else
    {
      error_flag = evolveFun(mem, final_time, y, &tret, normal);
    }

    states.emplace_back(n_vector_to_state(y, user_data.state_descs).deep_copy());
    outputFun();

    return RainshaftSolution(std::move(states), countFun());
  }

  static int postprocess_positive(sunrealtype, N_Vector y, void* user_data) {
    const RainshaftUserData& cast_data = *static_cast<RainshaftUserData *>(user_data);
    State state = n_vector_to_state(y, cast_data.state_descs);
    VarMut t = state.get_variable("T");
    VarMut q = state.get_variable("q");
    VarMut nr = state.get_variable("nr");
    VarMut qr = state.get_variable("qr");

    for (std::size_t i = 0; i != t.size(); ++i) {
      if (qr[i] < cast_data.constants.qsmall) {
        // Note that for the "original" P3 method applied to the rainshaft model,
        // which this class was designed for, the only way to get significant
        // negative qr is from excessive evaporation (a source of q), so this
        // should never drive q negative. Notably, this is only true so long as
        // the CFL condition is not violated in the sedimentation of qr (which
        // should not be a problem so long as sed_cfl_integrator is used).
        q[i] += qr[i];
        t[i] -= qr[i] * cast_data.constants.latvap / cast_data.constants.cp;
        qr[i] = 0.;
      }
      // nr is not conserved by local processes, but we do apply size limiters.
      nr[i] = cast_data.size_limiters.limited_nr(nr[i], qr[i]);
      // Note: This t limiter is not in the original P3 code, which may silently
      // accept negative t or crash on negative t, depending on treatment of
      // floating point exceptions and where the negative t occurs.
      t[i] = std::max(t[i], 0.);
      // In the rainshaft model, large negative q values should not be produced
      // since there is no "sink" except nudging (which should only be used with
      // time scales >> dt), but the below will violate energy/mass conservation
      // if more processes are added without adapting more of P3's limiter logic.
      q[i] = std::max(q[i], 0.);
    }

    return 0;
  }

  // Create an N_Vector from a StateConst. To preserve const correctness, this
  // makes a copy of the data.
  static N_Vector view_to_n_vector(const sundials::Context &sun_ctxt, const spaecies::VariableArrayView<const double> &view)
  {
    N_Vector y = N_VNew_Serial(view.size(), sun_ctxt);
    std::copy_n(view.data(), view.size(), N_VGetArrayPointer(y));
    return y;
  }

  N_Vector fill_abs_tol_vector(N_Vector abs_tol) const {
    double * const tol_data = N_VGetArrayPointer(abs_tol);
    const std::size_t nz = user_data.grid.nlev;
    // std::fill_n(tol_data, nz, 1.e-1);
    // std::fill_n(&tol_data[nz], nz, 1.e-3);
    // std::fill_n(&tol_data[2 * nz], nz, 1.e-4);
    // std::fill_n(&tol_data[3 * nz], nz, 1.e-12);
    // std::fill_n(tol_data, nz, 1.e-4);
    // std::fill_n(&tol_data[nz], nz, 1.e-6);
    // std::fill_n(&tol_data[2 * nz], nz, 1.e-7);
    // std::fill_n(&tol_data[3 * nz], nz, 1.e-15);
    // std::fill_n(tol_data, nz, 1.e-8);
    // std::fill_n(&tol_data[nz], nz, 1.e-10);
    // std::fill_n(&tol_data[2 * nz], nz, 1.e-11);
    // std::fill_n(&tol_data[3 * nz], nz, 1.e-19);
    std::fill_n(tol_data, nz, 0.e-1);
    std::fill_n(&tol_data[nz], nz, 0.e-3);
    std::fill_n(&tol_data[2 * nz], nz, 0.e-4);
    std::fill_n(&tol_data[3 * nz], nz, 0.e-12);
    return abs_tol;
  }

private:
  // The state and tendency conversions below are currently duplicated, but
  // in the future the arguments to the state and tendency constructors will be different...
  // Note that in both cases the N_Vector is treated as the owner of the data, and must
  // outlive the state/tendency objects (unless the deep_copy method is used to get a newly
  // allocated copy).

  static State n_vector_to_state(N_Vector y, const VarDescList& var_descs)
  {
    return {var_descs, N_VGetArrayPointer(y)};
  }

  static Tendency n_vector_to_tendency(N_Vector y, const VarDescList& var_descs)
  {
    return {var_descs, N_VGetArrayPointer(y)};
  }
};

#endif // SUNDIALS_INTEGRATOR_HPP
