#ifndef SPAECIES_STATE_HPP
#define SPAECIES_STATE_HPP

#include "variable_array_view.hpp"

namespace spaecies {

template <class T>
class State : public VariableArrayView<T> {
public:
  State(std::initializer_list<VarDescPtr> var_descs)
    : VariableArrayView<T>(var_descs) {
  }
  State(const std::vector<VarDescPtr>& var_descs)
    : VariableArrayView<T>(var_descs) {
  }
  // It is likely that some future functionality will be desired here, but for
  // now, we produce this specialization just to document the purpose of the
  // type and make sure that mixing up arrays with different purposes results in
  // a compile-time error.
  State(const std::vector<VarDescPtr>& var_descs, T* data_ptr)
    : VariableArrayView<T>(var_descs, data_ptr) {
  };
  State<T> deep_copy() const {
    State<T> copy(this->var_descs());
    this->copy_data_to_location(copy.data());
    return copy;
  }
};

}

#endif // SPAECIES_STATE_HPP
