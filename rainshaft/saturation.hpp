#ifndef SATURATION_HPP
#define SATURATION_HPP
#include "rainshaft_constants.hpp"
#include "derivatives.hpp"
#include <cmath>

// SPS: Should consider making this an abstract class with different formulas,
// since there are a number of parameterizations in use and they should be
// consistent across the model.

class SaturationFormulae {

public:

  // Initialize from constants.
  SaturationFormulae(const RainshaftConstants& constants);

  // Saturation vapor pressure over liquid at a given temperature.
  template <bool WithGrad = false>
  Val<WithGrad, 1> svp_liquid(const double temperature) const {
    const auto log_t = std::log(temperature);
    const auto recip_t = 1. / temperature;
    const auto tanh_fac = std::tanh(0.0415 * (temperature - 218.8));
    const auto scaled_recip_t_1 = 1331.22 * recip_t;
    const auto tanh_coeff = 53.878 - scaled_recip_t_1 - 9.44523 * log_t + 0.014025 * temperature;
    const auto tanh_term = tanh_coeff * tanh_fac;
    const auto scaled_recip_t_2 = 6763.22 * recip_t;
    const auto log_esl = 54.842763 - scaled_recip_t_2 - 4.21 * log_t
      + 0.000367*temperature + tanh_term;
    const auto svp = std::exp(log_esl);

    if constexpr (WithGrad) {
      const auto tanh_fac_grad = 0.0415 * tanh_coeff * (1.0 - pow(tanh_fac, 2)) + tanh_fac * (0.014025 + (scaled_recip_t_1 - 9.44523) * recip_t);
      return {svp, {
        svp * (0.000367 + (scaled_recip_t_2 - 4.21) * recip_t + tanh_fac_grad)
      }};
    } else {
      return svp;
    }
  }

  // Saturation specific humidity (over mass of dry air) as a function of
  // temperature and dry pressure.
  template <bool WithGrad = false>
  Val<WithGrad, 1> q_sat_dry(const double temperature, const double pressure_dry) const {
    const auto esl = svp_liquid<WithGrad>(temperature);
    const auto fac = epsilon_h2o / pressure_dry;
    const auto qds = fac * get_val(esl);

    if constexpr (WithGrad) {
      const auto [esl_grad] = get_grad(esl);
      return {qds, {fac * esl_grad}};
    } else {
      return qds;
    }
  }

private:

  // Copy of ratio of water vapor molecular mass to that of dry air.
  // Stored mainly as a convenience to not require a RainshaftConstants as an
  // argument to a bunch of functions.
  double epsilon_h2o;

};

#endif // SATURATION_HPP
