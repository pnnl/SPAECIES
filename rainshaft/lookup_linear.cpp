#include "lookup_linear.hpp"
#include <cmath>
#include <algorithm>

RealGrad<1> LookupLinear::lookup_value(double x) const {
  // Binary search to find first bound > x
  auto itr = std::upper_bound(range_bounds.cbegin(), range_bounds.cend(), x);

  // If x outside the range, return the endpoint value with 0 slope
  if (itr == range_bounds.cbegin()) {
    return {*itr, {0.0}};
  }
  if (itr == range_bounds.cend()) {
    return {*std::prev(itr), {0.0}};
  }

  // Index of largest range bound <= x
  const std::size_t range_idx = std::distance(range_bounds.cbegin(), std::prev(itr));
  const double spacing = spacings[range_idx];
  double int_part;
  const double frac_part = std::modf((x - *std::prev(itr)) / spacing, &int_part);
  const std::size_t cell_idx = static_cast<std::vector<double>::size_type>(int_part);
  const double y_left = y_table[range_idx][cell_idx];
  const double y_diff = y_table[range_idx][cell_idx + 1] - y_left;
  return {y_left + y_diff * frac_part, {y_diff / spacing}};
}
