#include "domain.hpp"

#include <algorithm>

namespace spaecies {

DimensionPtr Domain::add_dimension(const std::string& name, std::size_t size) {
  DimensionPtr new_dim = std::make_shared<Dimension>(name, size);
  dimensions.push_back(new_dim);
  return new_dim;
}

DimensionPtr Domain::get_dimension(const std::string& name) {
  auto matches_name = [&](DimensionPtr& dim_ptr) {
    return dim_ptr->name == name;
  };
  DimensionPtr dim_ptr = *std::find_if(dimensions.begin(), dimensions.end(),
                                       matches_name);
  return dim_ptr;
}

VarDescPtr Domain::add_var_desc(const std::string name,
                                VariableType type,
                                const std::vector<DimensionPtr> dimensions,
                                const std::string units,
                                VariableConstantStatus constant_status,
                                const std::optional<const std::string>& description,
                                const std::optional<const std::string>& standard_name) {
  VarDescPtr var_desc = std::make_shared<VariableDescriptor>(name, type, dimensions,
                                                             units, constant_status,
                                                             description, standard_name);
  var_descs.push_back(var_desc);
  return var_desc;
}

// Retrieve the dimension with a given size.
VarDescPtr Domain::get_var_desc(const std::string& name) {
  auto matches_name = [&](VarDescPtr& var_desc) {
    return var_desc->name == name;
  };
  VarDescPtr var_desc = *std::find_if(var_descs.begin(), var_descs.end(),
                                      matches_name);
  return var_desc;
}

}
