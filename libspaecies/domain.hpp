#ifndef SPAECIES_DOMAIN_HPP
#define SPAECIES_DOMAIN_HPP

#include <vector>

#include "dimension.hpp"

namespace spaecies {

class Domain {
public:
  // Add a dimension to this domain, returning a pointer to the constructed
  // dimension.
  DimensionPtr add_dimension(const std::string& name, std::size_t size);
  // Retrieve the dimension with a given size.
  DimensionPtr get_dimension(const std::string& name);
private:
  std::vector<DimensionPtr> dimensions;
};

}

#endif // SPAECIES_DOMAIN_HPP
