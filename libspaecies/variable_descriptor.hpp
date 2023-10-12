#ifndef SPAECIES_VARIABLE_DESCRIPTOR_HPP
#define SPAECIES_VARIABLE_DESCRIPTOR_HPP

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "dimension.hpp"

namespace spaecies
{

enum VariableType {
  Float64Type,
  Float32Type,
  Int64Type,
  Int32Type,
  BoolType
};

enum VariableConstantStatus {
  IsConstant,
  IsNotConstant
};

class VariableDescriptor {
public:

  VariableDescriptor(const std::string name,
                     VariableType type,
                     const std::vector<DimensionPtr> dimensions,
                     const std::string units,
                     VariableConstantStatus constant_status=IsNotConstant,
                     const std::optional<const std::string>& description=std::optional<const std::string>(),
                     const std::optional<const std::string>& standard_name=std::optional<const std::string>());
  // Short (but unique) name identifying this variable.
  const std::string name;
  // Type of the variable.
  const VariableType type;
  // List of dimensions defining the domain of this variable.
  const std::vector<DimensionPtr> dimensions;
  // Units associated with the variable as a string.
  const std::string units;
  // Whether the variable is constant over time.
  const VariableConstantStatus constant_status;
  // Metadata describing the meaning of this variable.
  const std::optional<const std::string> description;
  // "Standard" name that may be used by external tools.
  const std::optional<const std::string> standard_name;
  // Total number of values associated with the variable, i.e. the product of
  // all dimension sizes.
  std::size_t size();

};

}

#endif // SPAECIES_VARIABLE_DESCRIPTOR_HPP
