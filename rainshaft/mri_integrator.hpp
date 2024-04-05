#ifndef MRI_INTEGRATOR_HPP
#define MRI_INTEGRATOR_HPP
#include "sundials_integrator.hpp"

class MRIIntegrator : public SundialsIntegrator {

public:

  int fastOrder;
  int slowOrder;
  double dt;
  double dt_slowfac;

  MRIIntegrator(const RainshaftConstants* constants,
                     const int fast_order,
                     const int slow_order,
                     const double dt_in,
                     const double dt_slowfacin,
                     const RainshaftGrid* grid,
                     const RainshaftProcess* process_exp,
                     const RainshaftProcess* process_imp,
                     const RainshaftProcess* process_all,
                     sundials::Context *sun_ctxt);



  RainshaftSolution integrate(double initial_time,
                              double final_time,
                              const RainshaftState& initial_state) const;

};

#endif // MRI_INTEGRATOR_HPP
