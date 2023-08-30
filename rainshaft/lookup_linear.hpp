#ifndef LOOKUP_LINEAR_HPP
#define LOOKUP_LINEAR_HPP
#include <vector>

#include "lookup_table.hpp"

// Lookup table using linear interpolation.
// Assumes that the input ("x") range can be split into sub-ranges, each
// of which has a uniform spacing of x values within them.
// Extrapolated values use the nearest table value.
class LookupLinear : public LookupTable {

public:

  LookupLinear(std::vector<double> x_range_bounds,
               std::vector<double> x_spacings,
               std::vector<double> y_values);

  // Look up a table value.
  double lookup_value(double x) const;

private:

  // Bounds of x-value subranges.
  const std::vector<double> range_bounds;

  // Spacings between x-values in each subrange.
  const std::vector<double> spacings;

  // Index offsets for each subrange
  const std::vector<std::size_t> offsets;

  // Table of x-values to refer to for saving time.
  const std::vector<double> x_table;

  // Table of y-values to refer to.
  const std::vector<double> y_table;

};

#endif // LOOKUP_LINEAR_HPP
