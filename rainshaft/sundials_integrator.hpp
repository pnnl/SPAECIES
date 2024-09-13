#ifndef SUNDIALS_INTEGRATOR_HPP
#define SUNDIALS_INTEGRATOR_HPP
#include <array>
#include <stdexcept>
#include <limits>
#include "rainshaft_integrator.hpp"
#include "sundials/sundials_context.hpp"
#include "nvector/nvector_serial.h"

#include "rainshaft_integrator.hpp"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_process.hpp"
#include "rainshaft_solution.hpp"
#include "rainshaft_types.hpp"

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
    RainshaftUserData *cast_data = (RainshaftUserData *)user_data;
    const StateConst state = n_vector_to_state(y, cast_data->state_descs);
    RainshaftDerivedVars dvars = RainshaftDerivedVars(cast_data->constants,
                                                      cast_data->grid,
                                                      state);
    // Zero out ydot so that we don't have to remember to zero every value in calc_tend.
    N_VConst(0., ydot);
    const Tendency tend = n_vector_to_tendency(ydot, cast_data->tend_descs);
    cast_data->processes[PARTITION]->calc_tend(cast_data->constants,
                                               cast_data->grid,
                                               state, dvars, tend);
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
  decltype(&rainshaft_f<PARTITION>) create_f() const
  {
    if (user_data.processes[PARTITION] == nullptr)
    {
      return nullptr;
    }
    else
    {
      return rainshaft_f<PARTITION>;
    }
  }

  template <typename E, typename C, typename Mode>
  RainshaftSolution evolve(E evolveFun,
                           C countFun,
                           void *mem,
                           const double final_time,
                           const N_Vector y0,
                           const Mode normal,
                           const Mode one_step) const
  {
    std::vector<StateConst> states;

    sunrealtype tret = -std::numeric_limits<sunrealtype>::infinity();
    if (steps_per_output > 0)
    {
      states.emplace_back(n_vector_to_state(y0, user_data.state_descs));
      N_Vector y_out = N_VClone(y0);
      bool output_this_iter = true;
      for (int i = 0; tret < final_time; i++)
      {
        evolveFun(mem, final_time, y_out, &tret, one_step);
        output_this_iter = (i+1) % steps_per_output == 0;
        if (output_this_iter)
        {
          // Copy because (a) y_out is reused between iterations, and
          // (b) it owns its data.
          states.emplace_back(n_vector_to_state(y_out, user_data.state_descs).deep_copy());
        }
      }
      // Output last iteration if this was not already done.
      if (!output_this_iter)
      {
        states.emplace_back(n_vector_to_state(y_out, user_data.state_descs).deep_copy());
      }
      N_VDestroy(y_out);
    }
    else
    {
      // If we're only evolving and saving output once, place output directly
      // into the state vector rather than copying.
      State final_state(user_data.state_descs); // empty state for output
      N_Vector y_out = view_to_n_vector(sun_ctxt, final_state); // y_out does not own memory
      evolveFun(mem, final_time, y_out, &tret, normal);
      states.emplace_back(std::move(final_state)); // Transfer owning state to solution vector.
      N_VDestroy(y_out);
    }

    return RainshaftSolution(std::move(states), countFun());
  }

  // Note that the N_Vector does not take ownership of the data and will not free it.
  // The N_Vector must not be permitted to live longer than the view.
  static N_Vector view_to_n_vector(const sundials::Context &sun_ctxt, const spaecies::VariableArrayView<double> &view)
  {
    sunindextype num_variables = view.size();
    N_Vector y = N_VMake_Serial(num_variables, view.data(), sun_ctxt);
    return y;
  }

  // Note (as above) that the N_Vector will not own the data, but rather will have a
  // non-const pointer to data that should be treated as const.
  //
  // There's no way to force SUNDIALS to respect this const-ness, so this is only safe
  // for cases where we know SUNDIALS won't modify the variables, e.g. for the initial
  // conditions.
  static N_Vector state_to_y0(const sundials::Context &sun_ctxt, const StateConst &view)
  {
    sunindextype num_variables = view.size();
    N_Vector y = N_VMake_Serial(num_variables, const_cast<double*>(view.data()), sun_ctxt);
    return y;
  }

  // The state and tendency conversions below are currently duplicated, but
  // in the future the arguments to the state and tendency constructors will be different...
  // Note that in both cases the N_Vector is treated as the owner of the data, and must
  // outlive the state/tendency objects (unless the deep_copy method is used to get a newly
  // allocated copy).

  static State n_vector_to_state(N_Vector y, const VarDescList& var_descs)
  {
    sunrealtype* ydata = N_VGetArrayPointer(y);
    return State(var_descs, &ydata[0]);
  }

  static Tendency n_vector_to_tendency(N_Vector y, const VarDescList& var_descs)
  {
    sunrealtype* ydata = N_VGetArrayPointer(y);
    return Tendency(var_descs, &ydata[0]);
  }
};

#endif // SUNDIALS_INTEGRATOR_HPP
