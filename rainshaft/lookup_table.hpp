#ifndef LOOKUP_TABLE_HPP
#define LOOKUP_TABLE_HPP
#include <vector>
#include <cstddef>

// Abstract base class for lookup tables.
class LookupTable {

public:

  virtual double lookup_value(double x) const = 0;

  static std::vector<std::size_t> subrange_offsets(std::vector<double> x_range_bounds,
                                                   std::vector<double> x_spacings);

  static std::vector<double> calc_x_values(std::vector<double> x_range_bounds,
                                           std::vector<double> x_spacings);

};

#endif // LOOKUP_TABLE_HPP
