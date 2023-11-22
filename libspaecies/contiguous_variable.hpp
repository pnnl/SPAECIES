#ifndef SPAECIES_CONTIGUOUS_VARIABLE_HPP
#define SPAECIES_CONTIGUOUS_VARIABLE_HPP

#include "variable.hpp"

namespace spaecies {

template<class T>
class ContiguousVariable : public Variable<T> {
public:
  // Construct a ContiguousVariable from a description and pointer.
  ContiguousVariable(const VarDescPtr var_desc, T* data_ptr)
    : Variable<T>(var_desc), data_ptr(data_ptr) {
  }
  // SPS: Fix unsafe access issue here.
  // Construct a ContiguousVariable from a description and const pointer.
  ContiguousVariable(const VarDescPtr var_desc, const T* data_ptr)
    : Variable<T>(var_desc), data_ptr((T*) data_ptr) {
  }
  // Access a scalar value of the flattened variable.
  inline T& operator[](std::size_t idx) {
    return data_ptr[idx];
  }
  inline const T operator[](std::size_t idx) const {
    return data_ptr[idx];
  }

private:
  T *data_ptr;

};

}

#endif // SPAECIES_CONTIGUOUS_VARIABLE_HPP
