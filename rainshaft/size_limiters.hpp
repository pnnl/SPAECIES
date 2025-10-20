#ifndef SIZE_LIMITERS_HPP
#define SIZE_LIMITERS_HPP

#include "rainshaft_constants.hpp"

class SizeLimiters {

public:

  SizeLimiters(const RainshaftConstants& constants,
               double min_diameter, double max_diameter, double mu);

  double limited_nr(double nr, double qr) const;

private:

  const double min_nr_fac, max_nr_fac;

};

#endif // SIZE_LIMITERS_HPP
