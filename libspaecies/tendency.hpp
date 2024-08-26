#ifndef SPAECIES_TENDENCY_HPP
#define SPAECIES_TENDENCY_HPP

#include "variable_array_view.hpp"

namespace spaecies {

template <class T>
class Tendency : public VariableArrayView<T> {
public:
  Tendency(std::initializer_list<VarDescPtr> var_descs)
    : VariableArrayView<T>(var_descs) {
  }
  Tendency(const std::vector<VarDescPtr>& var_descs)
    : VariableArrayView<T>(var_descs) {
  }
  // It is likely that some future functionality will be desired here, but for
  // now, we produce this specialization just to document the purpose of the
  // type and make sure that mixing up arrays with different purposes results in
  // a compile-time error.
  Tendency(const std::vector<VarDescPtr>& var_descs, T* data_ptr)
    : VariableArrayView<T>(var_descs, data_ptr) {
  };
  static const Tendency<T> make_const(const std::vector<VarDescPtr>& var_descs, const T* data_ptr) {
    return Tendency(var_descs, const_cast<T*>(data_ptr));
  }
  Tendency<T> deep_copy() const {
    Tendency<T> copy(this->var_descs());
    this->copy_data_to_location(copy.data());
    return copy;
  }
};

}

#endif // SPAECIES_TENDENCY_HPP
