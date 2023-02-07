#ifndef RAINSHAFT_GRID_HPP
#define RAINSHAFT_GRID_HPP
#include <cstddef>
#include <vector>

#include "rainshaft_constants.hpp"

class RainshaftGrid {

public:

  // Constructor from interface values
  RainshaftGrid(const std::vector<double>& p_int_vec,
                const std::vector<double>& z_int_vec);

  // Number of cell interfaces
  const std::size_t nlev;

  // Pressure at cell interfaces (Pa)
  const std::vector<double> p_int;
  // Pressure at cell "midpoints" (Pa)
  const std::vector<double> p_mid;
  // Altitude at cell interfaces (m)
  const std::vector<double> z_int;
  // Cell heights (m)
  const std::vector<double> dz;

};

RainshaftGrid make_e3sm_like_grid(RainshaftConstants constants,
                                  double model_top, double srf_pres,
                                  double srf_temp, double lapse_rate);

#endif // RAINSHAFT_GRID_HPP
