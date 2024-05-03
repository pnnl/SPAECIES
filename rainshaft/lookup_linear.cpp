#include "lookup_linear.hpp"
#include <cstddef>
#include <stdexcept>

LookupLinear::LookupLinear(std::vector<double> x_range_bounds,
                           std::vector<double> x_spacings,
                           std::vector<double> y_values)
  : range_bounds(x_range_bounds), spacings(x_spacings),
    offsets(subrange_offsets(x_range_bounds, x_spacings)),
    x_table(calc_x_values(x_range_bounds, x_spacings)),
    y_table(y_values) {
  // SPS: Different prerequisites to check:
  //  - Size of x_range_bounds is one more than size of x_spacings.
  //  - x_range_bounds is strictly monotonically increasing.
  //  - Each difference between bounds is (within a narrow tolerance)
  //    a multiple of the corresponding spacing.
  //  - Number of y_values matches the number of x values.
}

double LookupLinear::lookup_value(double x) const {
  std::size_t num_subranges = spacings.size();
  if (x < range_bounds[0]) {
    return y_table[0];
  }
  for (std::size_t isr = 0; isr != num_subranges; ++isr) {
    if (x < range_bounds[isr+1]) {
      std::size_t low_ind = (x - range_bounds[isr]) / spacings[isr] + offsets[isr];
      double frac_part = (x - x_table[low_ind]) / (x_table[low_ind+1] - x_table[low_ind]);
      return y_table[low_ind]
        + frac_part * (y_table[low_ind+1] - y_table[low_ind]);
    }
  }
  return y_table[offsets[num_subranges]];
}
