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
    const PartitionArray processes;
    const VarDescList state_descs;
    const VarDescList tend_descs;
  };

  template <int PARTITION>
  static int rainshaft_f(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
  {
    auto *cast_data = static_cast<RainshaftUserData *>(user_data);
    const StateConst state = n_vector_to_state(y, cast_data->state_descs);
    const auto dvars = RainshaftDerivedVars(cast_data->constants,
                                                      cast_data->grid,
                                                      state);
    // Zero out ydot so that we don't have to remember to zero every value in calc_tend.
    N_VConst(0., ydot);
    Tendency tend = n_vector_to_tendency(ydot, cast_data->tend_descs);
    cast_data->processes[PARTITION]->calc_tend(cast_data->constants,
                                               cast_data->grid,
                                               state, dvars, tend);
    return 0;
  }

  template <int PARTITION>
  static int rainshaft_jac_prod(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp)
  {
    N_VConst(0, Jv);
    auto *cast_data = static_cast<RainshaftUserData *>(user_data);
    const StateConst state = n_vector_to_state(y, cast_data->state_descs);
    const auto dvars = RainshaftDerivedVars(cast_data->constants,
                                      cast_data->grid,
                                      state);
    cast_data->processes[PARTITION]->calc_tend_jac_prod(cast_data->constants,
                                                        cast_data->grid,
                                                        state, dvars,
                                                        N_VGetArrayPointer(v),
                                                        N_VGetArrayPointer(Jv));
    return 0;
  }

  template <int PARTITION>
  static int rainshaft_jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
  {
    SUNMatZero(Jac);
    auto *cast_data = static_cast<RainshaftUserData *>(user_data);
    const StateConst state = n_vector_to_state(y, cast_data->state_descs);
    const auto dvars = RainshaftDerivedVars(cast_data->constants,
                                      cast_data->grid,
                                      state);

    RainshaftProcess::Matrix mat = [Jac](const auto i, const auto j) -> auto & {
      return SM_ELEMENT_D(Jac, i, j);
    };

    cast_data->processes[PARTITION]->calc_tend_jac(cast_data->constants,
                                                   cast_data->grid,
                                                   state, dvars,
                                                   mat);
    return 0;
  }

  static void handle_error(int line,
                           const char *func,
                           const char *file,
                           const char *msg,
                           SUNErrCode err_code,
                           void *err_user_data,
                           SUNContext ctx)
  {
    throw std::logic_error(msg);
  }

public:
  SundialsIntegrator(const RainshaftConstants &constants,
                     const RainshaftGrid &grid,
                     const PartitionArray processes,
                     const VarDescList& state_descs,
                     const VarDescList& tend_descs,
                     const int steps_per_output)
      : user_data{constants, grid, processes, state_descs, tend_descs}, steps_per_output(steps_per_output)
  {
    SUNContext_PushErrHandler(sun_ctxt, SundialsIntegrator::handle_error, nullptr);
  }

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
  auto create_jac_prod() const
  {
    return (user_data.processes[PARTITION] == nullptr) ? nullptr : rainshaft_jac_prod<PARTITION>;
  }

  template <int PARTITION>
  auto create_jac() const
  {
    return (user_data.processes[PARTITION] == nullptr) ? nullptr : rainshaft_jac<PARTITION>;
  }

  template <typename E, typename C, typename O, typename Mode>
  RainshaftSolution evolve(E evolveFun,
                           C countFun,
                           O outputFun,
                           void *mem,
                           const double final_time,
                           const N_Vector y,
                           const Mode normal,
                           const Mode one_step) const
  {
    std::vector<StateConst> states;
    auto tret = -std::numeric_limits<sunrealtype>::infinity();
    if (steps_per_output > 0)
    {
      for (int i = 0; tret < final_time; i++)
      {
        if (i % steps_per_output == 0) {
          states.emplace_back(n_vector_to_state(y, user_data.state_descs).deep_copy());
          outputFun();
        }
        evolveFun(mem, final_time, y, &tret, one_step);
      }
    }
    else
    {
      evolveFun(mem, final_time, y, &tret, normal);
    }

    states.emplace_back(n_vector_to_state(y, user_data.state_descs).deep_copy());
    outputFun();

    return RainshaftSolution(std::move(states), countFun());
  }

  // Create an N_Vector from a StateConst. To preserve const correctness, this
  // makes a copy of the data.
  static N_Vector state_to_n_vector(const sundials::Context &sun_ctxt, const StateConst &view)
  {
    N_Vector y = N_VNew_Serial(view.size(), sun_ctxt);
    std::copy_n(view.data(), view.size(), N_VGetArrayPointer(y));
    return y;
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
