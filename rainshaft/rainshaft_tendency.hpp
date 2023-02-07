#ifndef RAINSHAFT_TENDENCY_HPP
#define RAINSHAFT_TENDENCY_HPP
#include <vector>

class RainshaftTendency {

public:

  // Constructor from initial t, q, nr, and qr rates of change.
  RainshaftTendency(const std::vector<double>& t_vec,
                    const std::vector<double>& q_vec,
                    const std::vector<double>& nr_vec,
                    const std::vector<double>& qr_vec);

  // Temperature rate (K/s)
  const std::vector<double> t_tend;
  // Water vapor mass mixing ratio rate (kg/kg/s)
  const std::vector<double> q_tend;
  // Rain number concentration rate (#/kg/s)
  const std::vector<double> nr_tend;
  // Rain mass mixing ratio rate (kg/kg/s)
  const std::vector<double> qr_tend;
};

#endif // RAINSHAFT_TENDENCY_HPP
