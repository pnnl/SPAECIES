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
  // Cast to const
  operator Tendency<const T>() {
    return Tendency<const T>((VariableArrayView<T>)*this);
  }
  // Allow above use of private constructor to work.
  friend class Tendency<std::remove_const_t<T>>;
  // It is likely that some future functionality will be desired here, but for
  // now, we produce this specialization just to document the purpose of the
  // type and make sure that mixing up arrays with different purposes results in
  // a compile-time error.
  Tendency(const std::vector<VarDescPtr>& var_descs, T* data_ptr)
    : VariableArrayView<T>(var_descs, data_ptr) {
  };
  Tendency<T> deep_copy() const {
    Tendency<T> copy(this->var_descs());
    this->copy_data_to_location(copy.data());
    return copy;
  }
private:
  // Making this private because we don't want to guarantee that a Tendency
  // can always be constructed from a VariableArrayView alone.
  Tendency(const VariableArrayView<T>& view) : VariableArrayView<T>(view) {
  }
};

}

#endif // SPAECIES_TENDENCY_HPP
