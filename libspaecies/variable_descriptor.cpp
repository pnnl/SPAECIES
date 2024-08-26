#include "variable_descriptor.hpp"

#include <functional>
#include <numeric>

#include "exceptions.hpp"

namespace spaecies
{

std::string spaecies_type_name(VariableType type) {
  switch (type) {
  case Float64Type:
    return FLOAT64NAME;
  case Float32Type:
    return FLOAT32NAME;
  case Int64Type:
    return INT64NAME;
  case Int32Type:
    return INT32NAME;
  case BoolType:
    return BOOLNAME;
  case InvalidType:
    return INVALIDNAME;
  }
  throw UnreachableException("invalid VariableType in spaecies_type_name");
}

VariableDescriptor::VariableDescriptor(const std::string& name,
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
