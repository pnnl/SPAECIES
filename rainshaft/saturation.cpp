#include "saturation.hpp"

SaturationFormulae::SaturationFormulae(const RainshaftConstants& constants)
  : epsilon_h2o(constants.epsilon_h2o) {
}

double SaturationFormulae::svp_liquid(double temperature) {
  double log_t = std::log(temperature);
  double recip_t = 1./temperature;
  double tanh_fac = std::tanh(0.0415 * (temperature - 218.8));
  double tanh_term = (53.878 - 1331.22 * recip_t - 9.44523 * log_t
                      + 0.014025 * temperature) * tanh_fac;
  double log_esl = 54.842763 - 6763.22*recip_t - 4.21 * log_t
    + 0.000367*temperature + tanh_term;
  return std::exp(log_esl);
}

double SaturationFormulae::wv_pressure_to_q_dry(double pressure_wv,
                                                double pressure_dry) {
  return epsilon_h2o * pressure_wv / pressure_dry;
}

double SaturationFormulae::q_sat_dry(double temperature, double pressure_dry) {
  double esl = svp_liquid(temperature);
  return wv_pressure_to_q_dry(esl, pressure_dry);
}
