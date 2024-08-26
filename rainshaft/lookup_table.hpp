#ifndef LOOKUP_TABLE_HPP
#define LOOKUP_TABLE_HPP
#include <vector>
#include <cstddef>
#include "derivatives.hpp"

// Abstract base class for lookup tables.
class LookupTable {

public:
  virtual ValGrad<1> lookup_value(double x) const = 0;
};

#endif // LOOKUP_TABLE_HPP
