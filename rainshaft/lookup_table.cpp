#include "lookup_table.hpp"
#include <cmath>

std::vector<std::size_t> LookupTable::subrange_offsets(std::vector<double> x_range_bounds,
                                                       std::vector<double> x_spacings) {
  std::size_t num_subranges = x_spacings.size();
  std::vector<std::size_t> offsets(num_subranges+1, 0.);
  for (std::size_t isr = 0; isr != num_subranges; ++isr) {
    double subrange_size = x_range_bounds[isr+1] - x_range_bounds[isr];
    std::size_t num_subrange_values = std::round(subrange_size / x_spacings[isr]);
    offsets[isr+1] = offsets[isr] + num_subrange_values;
  }
  return offsets;
}

std::vector<double> LookupTable::calc_x_values(std::vector<double> x_range_bounds,
                                               std::vector<double> x_spacings) {
  std::size_t num_subranges = x_spacings.size();
  std::vector<double> x_values;
  x_values.push_back(x_range_bounds[0]);
  for (std::size_t isr = 0; isr != num_subranges; ++isr) {
    double subrange_size = x_range_bounds[isr+1] - x_range_bounds[isr];
    std::size_t num_subrange_values = std::round(subrange_size / x_spacings[isr]);
    for (std::size_t isv = 0; isv != num_subrange_values; ++isv) {
      double ratio = ((double) (isv + 1)) / num_subrange_values;
      x_values.push_back(x_range_bounds[isr] + subrange_size * ratio);
    }
  }
  return x_values;
}
