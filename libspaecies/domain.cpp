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
  // SPS: need to cover case where variable is not found.
  return dim_ptr;
}

std::vector<DimensionPtr> Domain::get_dimensions(const std::vector<std::string>& names) {
  std::vector<DimensionPtr> dim_ptrs(names.size());
  std::transform(names.cbegin(), names.cend(),
                 dim_ptrs.begin(),
                 [&](const std::string name) { return get_dimension(name); });
  return dim_ptrs;
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
  // SPS: need to cover case where variable is not found.
  return var_desc;
}

std::vector<VarDescPtr> Domain::get_var_descs(const std::vector<std::string>& names) {
  std::vector<VarDescPtr> var_desc_ptrs(names.size());
  std::transform(names.cbegin(), names.cend(),
                 var_desc_ptrs.begin(),
                 [&](const std::string name) { return get_var_desc(name); });
  return var_desc_ptrs;
}

}
