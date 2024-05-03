#ifndef SUNDIALS_INTEGRATOR_HPP
#define SUNDIALS_INTEGRATOR_HPP
#include <array>
#include <stdexcept>
#include <limits>
#include "rainshaft_integrator.hpp"
#include "sundials/sundials_context.hpp"
#include "nvector/nvector_serial.h"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_state.hpp"
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
  };

  template <int PARTITION>
  static int rainshaft_f(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
  {
    // SPS: Should stop using std::vector to reduce copies and allocations.
    RainshaftState state = n_vector_to_state(y);
    RainshaftUserData *cast_data = (RainshaftUserData *)user_data;
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
                     const int steps_per_output)
      : user_data{constants, grid, processes}, steps_per_output(steps_per_output)
  {
    SUNContext_PushErrHandler(sun_ctxt, SundialsIntegrator::handle_error, nullptr);
  }

protected:
  const int steps_per_output;
  const sundials::Context sun_ctxt;
  const RainshaftUserData user_data;

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
                           const N_Vector y,
                           const Mode normal,
                           const Mode one_step) const
  {
    std::vector<RainshaftState> states;
    std::vector<RainshaftDerivedVars> dvars;

    sunrealtype tret = -std::numeric_limits<sunrealtype>::infinity();
    if (steps_per_output > 0)
    {
      for (int i = 0; tret < final_time; i++)
      {
        if (i % steps_per_output == 0)
        {
          const auto new_state = n_vector_to_state(y);
          dvars.emplace_back(user_data.constants, user_data.grid, new_state);
          states.push_back(new_state);
        }
        evolveFun(mem, final_time, y, &tret, one_step);
      }
    }
    else
    {
      evolveFun(mem, final_time, y, &tret, normal);
    }

    const auto new_state = n_vector_to_state(y);
    dvars.emplace_back(user_data.constants, user_data.grid, new_state);
    states.push_back(new_state);

    return RainshaftSolution(states, dvars, countFun());
  }

  static N_Vector state_to_n_vector(const sundials::Context &sun_ctxt, const RainshaftState &state)
  {
    sunindextype nz = state.t.size();
    // SPS: It should not assume 4 variables; move some of this to RainshaftState itself?
    sunindextype num_variables = nz * 4;
    N_Vector y = N_VNew_Serial(num_variables, sun_ctxt);
    sunrealtype *ydata = N_VGetArrayPointer_Serial(y);
    for (sunindextype j = 0; j != nz; ++j)
    {
      ydata[j] = state.t[j];
    }
    for (sunindextype j = 0; j != nz; ++j)
    {
      ydata[nz + j] = state.q[j];
    }
    for (sunindextype j = 0; j != nz; ++j)
    {
      ydata[2 * nz + j] = state.nr[j];
    }
    for (sunindextype j = 0; j != nz; ++j)
    {
      ydata[3 * nz + j] = state.qr[j];
    }
    return y;
  }

  // SPS: should reuse most of the equivalent state function here
  static void tend_to_n_vector(const RainshaftTendency &tend, N_Vector ydot)
  {
    sunindextype nz = tend.t_tend.size();
    sunrealtype *ydata = N_VGetArrayPointer_Serial(ydot);
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

  static RainshaftState n_vector_to_state(N_Vector y)
  {
    sunindextype nz = N_VGetLength(y) / 4;
    sunrealtype *ydata = N_VGetArrayPointer(y);
    std::vector<double> t(nz), q(nz), nr(nz), qr(nz);
    for (sunindextype j = 0; j != nz; ++j)
    {
      t[j] = ydata[j];
    }
    for (sunindextype j = 0; j != nz; ++j)
    {
      q[j] = ydata[nz + j];
    }
    for (sunindextype j = 0; j != nz; ++j)
    {
      nr[j] = ydata[2 * nz + j];
    }
    for (sunindextype j = 0; j != nz; ++j)
    {
      qr[j] = ydata[3 * nz + j];
    }
    return RainshaftState(t, q, nr, qr);
  }
};

#endif // SUNDIALS_INTEGRATOR_HPP
