#include "size_limiters.hpp"
#include <algorithm>

double diameter_to_nr_fac(const RainshaftConstants& constants, double diameter) {
  // Convert a diameter to the number which, when multiplied by qr, yields the
  // corresponding nr value for a gamma distribution with that mean diameter and
  // shape parameter mu.
  double lambdar = (constants.mur + 1.) / diameter;
  double mu_polynomial = (constants.mur + 1.) * (constants.mur + 2.) * (constants.mur + 3.);
  return lambdar*lambdar*lambdar*6. / (constants.pi * constants.rhow * mu_polynomial);
}

// Note that the *maximum* diameter sets the *minimum* nr limiter, and vice versa.
SizeLimiters::SizeLimiters(const RainshaftConstants& constants,
                           double min_diameter, double max_diameter)
  : min_nr_fac(diameter_to_nr_fac(constants, max_diameter)),
    max_nr_fac(diameter_to_nr_fac(constants, min_diameter))
{
}

double SizeLimiters::limited_nr(double nr, double qr) const {
  return std::max(std::min(nr, max_nr_fac * qr), min_nr_fac * qr);
}
