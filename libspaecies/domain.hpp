#ifndef SPAECIES_DOMAIN_HPP
#define SPAECIES_DOMAIN_HPP

#include <vector>

#include "dimension.hpp"
#include "variable_descriptor.hpp"

namespace spaecies {

class Domain {
public:
  // Add a dimension to this domain, returning a pointer to the constructed
  // dimension.
  DimensionPtr add_dimension(const std::string& name, std::size_t size);
  // Retrieve the dimension with a given size.
  DimensionPtr get_dimension(const std::string& name);
  // Add a variable to this domain, returning a pointer to the constructed
  // variable descriptor.
  VarDescPtr add_var_desc(const std::string name,
                          VariableType type,
                          const std::vector<DimensionPtr> dimensions,
                          const std::string units,
                          VariableConstantStatus constant_status=IsNotConstant,
                          const std::optional<const std::string>& description=std::optional<const std::string>(),
                          const std::optional<const std::string>& standard_name=std::optional<const std::string>());
  // Retrieve the dimension with a given size.
  VarDescPtr get_var_desc(const std::string& name);

private:
  // Array of all dimensions in the domain.
  std::vector<DimensionPtr> dimensions;
  // Array of all variable descriptors for variables in the domain.
  std::vector<VarDescPtr> var_descs;
};

}

#endif // SPAECIES_DOMAIN_HPP
