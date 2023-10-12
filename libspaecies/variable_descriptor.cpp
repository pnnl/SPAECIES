#include "variable_descriptor.hpp"

#include <functional>
#include <numeric>

namespace spaecies
{

VariableDescriptor::VariableDescriptor(const std::string name,
                                       VariableType type,
                                       const std::vector<DimensionPtr> dimensions,
                                       const std::string units,
                                       VariableConstantStatus constant_status,
                                       const std::optional<const std::string>& description,
                                       const std::optional<const std::string>& standard_name)
: name(name), type(type), dimensions(dimensions), units(units),
  constant_status(constant_status), description(description), standard_name(standard_name) {
}

std::size_t VariableDescriptor::size() {
  // Folding not available until C++23, so have to write loop explicitly.
  std::size_t var_size = 1;
  for (DimensionPtr d : dimensions) {
    var_size *= d->size;
  }
  return var_size;
}

}
