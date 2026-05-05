#include "variable_descriptor.hpp"

#include <functional>
#include <numeric>
#include <utility>

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

VariableDescriptor::VariableDescriptor(std::string name,
                                       VariableType type,
                                       std::vector<DimensionPtr> dimensions,
                                       std::string units,
                                       VariableUsage usage,
                                       std::optional<std::string> description,
                                       std::optional<std::string> standard_name)
: name(std::move(name)), type(type), dimensions(std::move(dimensions)), units(std::move(units)),
  usage(usage), description(std::move(description)), standard_name(std::move(standard_name)) {
}

std::size_t VariableDescriptor::size() {
  return std::accumulate(dimensions.begin(), dimensions.end(), 1,
                         [](std::size_t acc, const DimensionPtr& dim) { return acc * dim->size; });
}

}
