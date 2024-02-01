#ifndef RAINSHAFT_TENDENCY_JAC_HPP
#define RAINSHAFT_TENDENCY_JAC_HPP
#include <vector>

class RainshaftTendencyJac {

public:

  // Constructor from initial t, q, nr, and qr rates of change.
  RainshaftTendencyJac(const double* t_jac,
                    const double* q_jac,
                    const double* nr_jac,
                    const double* qr_jac);

  // Temperature rate (K/s)
  const double* t_tend_jac;
  // Water vapor mass mixing ratio rate (kg/kg/s)
  const double* q_tend_jac;
  // Rain number concentration rate (#/kg/s)
  const double* nr_tend_jac;
  // Rain mass mixing ratio rate (kg/kg/s)
  const double* qr_tend_jac;
};

#endif // RAINSHAFT_TENDENCY_JAC_HPP
