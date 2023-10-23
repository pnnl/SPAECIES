#ifndef SPAECIES_VARIABLE_HPP
#define SPAECIES_VARIABLE_HPP

#include "variable_descriptor.hpp"

namespace spaecies {

template<class T>
class Variable {
public:
  const VarDescPtr var_desc;
  // Access a scalar value of the flattened variable.
  virtual T& operator[](std::size_t idx) = 0;
  virtual const T operator[](std::size_t idx) const = 0;
  // Get the variable size.
  std::size_t size() {
    return var_desc->size();
  }

protected:
  Variable(const VarDescPtr var_desc) : var_desc(var_desc) {
  };

};

}

#endif // SPAECIES_VARIABLE_HPP
