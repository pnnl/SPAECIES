#ifndef SPAECIES_VARIABLE_ARRAY_HPP
#define SPAECIES_VARIABLE_ARRAY_HPP

#include <cstring>
#include <initializer_list>
#include <tuple>

#include "contiguous_variable_view.hpp"
#include "exceptions.hpp"
#include "variable_array_view.hpp"

namespace spaecies {

template <class T>
class VariableArray : public VariableArrayView<T> {
public:
  VariableArray(std::initializer_list<VarDescPtr> var_descs)
    : VariableArray(std::vector<VarDescPtr>(var_descs)) {
  }
  VariableArray(const std::vector<VarDescPtr>& var_descs)
    : VariableArrayView<T>(var_descs, new T[VariableArrayView<T>::calc_size(var_descs)]) {
  };
  VariableArray(const VariableArrayView<T>& view)
    : VariableArray(view.var_descs) {
    std::memcpy(this->data_ptr, view.data(), this->size * sizeof(T));
  }
  // Override copy constructor to prevent overlapping memory between arrays.
  VariableArray(const VariableArray<T>& other)
    : VariableArray((const VariableArrayView<T>&) other) {
  }
  // Move constructor.
  VariableArray(VariableArray<T>&& other)
    : VariableArrayView<T>(std::move(other)){
    other.data_ptr = nullptr;
  }
  ~VariableArray() {
    if (this->data_ptr) {
      delete[] this->data_ptr;
    }
  }
};

}

#endif // SPAECIES_VARIABLE_ARRAY_HPP
