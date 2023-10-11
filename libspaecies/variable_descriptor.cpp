#include "variable_descriptor.hpp"

namespace spaecies
{

VariableDescriptor::VariableDescriptor(const std::string name,
                                       VariableType type,
                                       const std::vector<DimensionPtr> dimensions,
                                       const std::string units)
  : name(name), type(type), dimensions(dimensions), units(units) {
}

std::size_t VariableDescriptor::size() {
  return dimensions[0]->size;
}

}
