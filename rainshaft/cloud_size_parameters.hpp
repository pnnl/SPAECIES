#ifndef CLOUD_SIZE_PARAMETERS_HPP
#define CLOUD_SIZE_PARAMETERS_HPP
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_types.hpp"

double calc_muc(double nc, double rho_dry);

double calc_lambdac(const RainshaftConstants& constants, double nc, double qc, double muc);

void limit_nc(const RainshaftConstants& constants,
              const RainshaftGrid& grid,
              State& state);

#endif // CLOUD_SIZE_PARAMETERS_HPP
