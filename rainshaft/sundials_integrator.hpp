#ifndef SUNDIALS_INTEGRATOR_HPP
#define SUNDIALS_INTEGRATOR_HPP
#include <array>
#include <stdexcept>
#include <limits>
#include "rainshaft_integrator.hpp"
#include "sundials/sundials_context.hpp"
#include "nvector/nvector_serial.h"

#include "spaecies.hpp"

#include "rainshaft_integrator.hpp"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_process.hpp"
#include "rainshaft_solution.hpp"
#include "rainshaft_tendency.hpp"

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
    const std::vector<spaecies::VarDescPtr> var_descs;
  };

  template <int PARTITION>
  static int rainshaft_f(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
  {
    RainshaftUserData *cast_data = (RainshaftUserData *)user_data;
    spaecies::VariableArrayView<double> state = n_vector_to_state(y, cast_data->var_descs);
    RainshaftDerivedVars dvars = RainshaftDerivedVars(cast_data->constants,
                                                      cast_data->grid,
                                                      state);
    RainshaftTendency tend = cast_data->processes[PARTITION]->calc_tend(cast_data->constants,
                                                                        cast_data->grid,
                                                                        state, dvars);
    tend_to_n_vector(tend, ydot);
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
                     const std::vector<spaecies::VarDescPtr>& var_descs,
                     const int steps_per_output)
      : user_data{constants, grid, processes, var_descs}, steps_per_output(steps_per_output)
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
                           const N_Vector y_init,
                           const Mode normal,
                           const Mode one_step) const
  {
    std::vector<spaecies::VariableArray<double>> states;
    std::vector<RainshaftDerivedVars> dvars;

    sunrealtype tret = -std::numeric_limits<sunrealtype>::infinity();
    if (steps_per_output > 0)
    {
      states.emplace_back(n_vector_to_state(y_init, user_data.var_descs));
      dvars.emplace_back(user_data.constants, user_data.grid, states.back());
      N_Vector y_out = N_VClone(y_init);
      bool output_this_iter = true;
      for (int i = 0; tret < final_time; i++)
      {
        evolveFun(mem, final_time, y_out, &tret, one_step);
        output_this_iter = (i+1) % steps_per_output == 0;
        if (output_this_iter)
        {
          states.emplace_back(n_vector_to_state(y_out, user_data.var_descs));
          dvars.emplace_back(user_data.constants, user_data.grid, states.back());
        }
      }
      // Output only if there isn't already an output for the last iteration.
      if (!output_this_iter)
      {
        states.emplace_back(n_vector_to_state(y_out, user_data.var_descs));
        dvars.emplace_back(user_data.constants, user_data.grid, states.back());
      }
      N_VDestroy(y_out);
    }
    else
    {
      // If we're only evolving and saving output once, place output directly
      // into the state vector rather than copying.
      states.emplace_back(user_data.var_descs); // empty state for output
      N_Vector y_out = state_to_n_vector(sun_ctxt, states.back());
      evolveFun(mem, final_time, y_out, &tret, normal);
      dvars.emplace_back(user_data.constants, user_data.grid, states.back());
      N_VDestroy(y_out);
    }

    return RainshaftSolution(std::move(states), std::move(dvars), countFun());
  }

  static N_Vector state_to_n_vector(const sundials::Context &sun_ctxt, const spaecies::VariableArrayView<double> &state)
  {
    sunindextype num_variables = state.size;
    // Unsafe cast is necessary because SUNDIALS does not have separate const
    // types for input vectors it will not change; we have to trust that the
    // caller and SUNDIALS will not change the result here.
    N_Vector y = N_VMake_Serial(num_variables, (double*) state.data(), sun_ctxt);
    return y;
  }

  // SPS: should reuse most of the equivalent state function here
  static void tend_to_n_vector(const RainshaftTendency &tend, N_Vector ydot)
  {
    sunindextype nz = tend.t_tend.size();
    auto ydata = N_VGetArrayPointer_Serial(ydot);
    for (sunindextype j = 0; j != nz; ++j)
    {
      ydata[j] = tend.t_tend[j];
    }
    for (sunindextype j = 0; j != nz; ++j)
    {
      ydata[nz + j] = tend.q_tend[j];
    }
    for (sunindextype j = 0; j != nz; ++j)
    {
      ydata[2 * nz + j] = tend.nr_tend[j];
    }
    for (sunindextype j = 0; j != nz; ++j)
    {
      ydata[3 * nz + j] = tend.qr_tend[j];
    }
  }

  static const spaecies::VariableArrayView<double> n_vector_to_state(N_Vector y, const std::vector<spaecies::VarDescPtr>& var_descs)
  {
    const sunrealtype* ydata = N_VGetArrayPointer(y);
    return spaecies::make_variable_array_view(var_descs, &ydata[0]);
  }
};

#endif // SUNDIALS_INTEGRATOR_HPP
