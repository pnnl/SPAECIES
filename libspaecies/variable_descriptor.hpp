#ifndef SPAECIES_VARIABLE_DESCRIPTOR_HPP
#define SPAECIES_VARIABLE_DESCRIPTOR_HPP

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "dimension.hpp"

namespace spaecies
{

// Variable types declared in descriptor.
enum VariableType {
  Float64Type,
  Float32Type,
  Int64Type,
  Int32Type,
  BoolType,
  InvalidType
};

// Simple lookup to connect C++ types to the dynamic descriptor types.
template <class T>
constexpr VariableType SPAECIES_TYPE = InvalidType;
template <>
inline constexpr VariableType SPAECIES_TYPE<double> = Float64Type;
template <>
inline constexpr VariableType SPAECIES_TYPE<float> = Float32Type;
template <>
inline constexpr VariableType SPAECIES_TYPE<int_least64_t> = Int64Type;
template <>
inline constexpr VariableType SPAECIES_TYPE<int_least32_t> = Int32Type;
template <>
inline constexpr VariableType SPAECIES_TYPE<bool> = BoolType;

const std::string FLOAT64NAME = "64-bit float";
const std::string FLOAT32NAME = "32-bit float";
const std::string INT64NAME = "64-bit integer";
const std::string INT32NAME = "32-bit integer";
const std::string BOOLNAME = "boolean";
const std::string INVALIDNAME = "[invalid type]";

std::string spaecies_type_name(VariableType type);

// Is the variable a constant in time?
enum VariableConstantStatus {
  IsConstant,
  IsNotConstant
};

class VariableDescriptor {
public:

  VariableDescriptor(const std::string& name,
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

// This typedef is convenient due to almost always using pointers to work with
// variable description objects.
typedef std::shared_ptr<VariableDescriptor> VarDescPtr;

}

#endif // SPAECIES_VARIABLE_DESCRIPTOR_HPP
