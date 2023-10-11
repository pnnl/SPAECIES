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

}
