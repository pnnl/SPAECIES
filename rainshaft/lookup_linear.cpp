#include "lookup_linear.hpp"
#include <cmath>
#include <algorithm>

ValGrad<1> LookupLinear::lookup_value(double x) const {
  auto itr = std::lower_bound(range_bounds.cbegin(), range_bounds.cend(), x);
  if (itr == range_bounds.begin()) {
    return {y_table.front().front(), {0.0}};
  }
  if (itr == range_bounds.end()) {
    return {y_table.back().back(), {0.0}};
  }

  const auto range_idx = std::distance(spacings.cbegin(), std::prev(itr));
  const auto spacing = spacings[range_idx];
  double int_part;
  const auto frac_part = std::modf((x - *std::prev(itr)) / spacing, &int_part);
  const auto cell_idx = static_cast<std::vector<double>::size_type>(int_part);
  const auto y_l = y_table[range_idx][cell_idx];
  const auto y_r = y_table[range_idx][cell_idx + 1];
  const auto slope = (y_r - y_l) / spacing;
  return {y_l + slope * frac_part, {slope}};
}
