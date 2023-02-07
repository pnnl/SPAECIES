#ifndef RAINSHAFT_STATE_HPP
#define RAINSHAFT_STATE_HPP
#include <vector>

class RainshaftState {

public:

  // Constructor from initial t, q, nr, and qr.
  RainshaftState(const std::vector<double>& t_vec,
                 const std::vector<double>& q_vec,
                 const std::vector<double>& nr_vec,
                 const std::vector<double>& qr_vec);

  // Temperature (K)
  const std::vector<double> t;
  // Water vapor mass mixing ratio (kg/kg)
  const std::vector<double> q;
  // Rain number concentration (#/kg)
  const std::vector<double> nr;
  // Rain mass mixing ratio (kg/kg)
  const std::vector<double> qr;
};

#endif // RAINSHAFT_STATE_HPP
