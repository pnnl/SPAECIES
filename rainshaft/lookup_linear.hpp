#ifndef LOOKUP_LINEAR_HPP
#define LOOKUP_LINEAR_HPP
#include <vector>
#include <algorithm>

#include "lookup_table.hpp"

// Lookup table using linear interpolation.
// Assumes that the input ("x") range can be split into sub-ranges, each
// of which has a uniform spacing of x values within them.
// Extrapolated values use the nearest table value.
class LookupLinear : public LookupTable {

public:
  template<typename F>
  LookupLinear(const std::vector<double> &range_bounds,
               const std::vector<double> &spacings,
               const F f) : range_bounds(range_bounds), spacings(spacings), y_table(create_table(range_bounds, spacings, f))
  {}

  // Look up a table value.
  RealGrad<1> lookup_value(double x) const;

private:
  using Table = std::vector<std::vector<double>>;

  template<typename F>
  static Table create_table(const std::vector<double> &range_bounds,
               const std::vector<double> &spacings,
               const F f)
  {
    Table tab(spacings.size());
    for (Table::size_type i = 0; i < spacings.size(); i++) {
      for (Table::value_type::size_type j = 0;; j++) {
        const double x = range_bounds[i] + j * spacings[i];
        tab[i].push_back(f(x));
        if (x >= range_bounds[i + 1]) {
          break;
        }
      }
    }

    return tab;
  }

  // Bounds of x-value subranges.
  const std::vector<double> range_bounds;

  // Spacings between x-values in each subrange.
  const std::vector<double> spacings;

  // Table of y-values to refer to.
  const Table y_table;

};

#endif // LOOKUP_LINEAR_HPP
