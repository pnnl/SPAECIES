#ifndef SPAECIES_CONTIGUOUS_VARIABLE_HPP
#define SPAECIES_CONTIGUOUS_VARIABLE_HPP

#include <type_traits>

#include "variable_view.hpp"

namespace spaecies {

template<class T> class ContiguousVariable;

template<class T>
ContiguousVariable<T> make_contiguous_variable(const VarDescPtr var_desc, T* data_ptr) {
  return ContiguousVariable(var_desc, data_ptr);
}

template<class T>
const ContiguousVariable<T> make_contiguous_variable(const VarDescPtr var_desc, const T* data_ptr) {
  return ContiguousVariable(var_desc, (T*) data_ptr);
}

template<class T>
class ContiguousVariable : public VariableView<T> {
public:
  // Construct a ContiguousVariable from a description and pointer.
  ContiguousVariable(const VarDescPtr var_desc, T* data_ptr)
    : VariableView<T>(var_desc), data_ptr(data_ptr) {
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
