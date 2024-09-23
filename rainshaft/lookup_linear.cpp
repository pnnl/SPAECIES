#include "lookup_linear.hpp"
#include <cmath>
#include <algorithm>

RealGrad<1> LookupLinear::lookup_value(double x) const {
  // Binary search to find first bound > x
  auto itr = std::upper_bound(range_bounds.cbegin(), range_bounds.cend(), x);

  // If x outside the range, return the endpoint value with 0 slope
  if (itr == range_bounds.cbegin()) {
    return {y_table.front().front(), {0.0}};
  }
  if (itr == range_bounds.cend()) {
    return {y_table.back().back(), {0.0}};
  }

  // Index of largest range bound <= x
  const auto range_idx = std::distance(range_bounds.cbegin(), std::prev(itr));
  const auto spacing = spacings[range_idx];
  double int_part;
  const auto frac_part = std::modf((x - *std::prev(itr)) / spacing, &int_part);
  const auto cell_idx = static_cast<std::vector<double>::size_type>(int_part);
  const auto y_left = y_table[range_idx][cell_idx];
  const auto y_diff = y_table[range_idx][cell_idx + 1] - y_left;
  return {y_left + y_diff * frac_part, {y_diff / spacing}};
}
