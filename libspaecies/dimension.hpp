#ifndef SPAECIES_DIMENSION_HPP
#define SPAECIES_DIMENSION_HPP

#include <cstddef>
#include <memory>
#include <string>

namespace spaecies
{

class Dimension {
public:

  Dimension(const std::string& name, std::size_t size);
  // Name of this dimension.
  const std::string name;
  // Size of this (discretized) dimension.
  const std::size_t size;

};

// This typedef is convenient due to almost always using pointers to work with
// dimension objects.
typedef std::shared_ptr<Dimension> DimensionPtr;

}

#endif // SPAECIES_DIMENSION_HPP
