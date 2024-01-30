#ifndef SPAECIES_VARIABLE_ARRAY_HPP
#define SPAECIES_VARIABLE_ARRAY_HPP

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
    : VariableArrayView<T>(var_descs, nullptr),
      data_array(this->size) {
    this->data_ptr = data_array.data();
  };
  const ContiguousVariableView<T> get_variable(const std::string& name) const {
    auto [var_desc, idx] = VariableArrayView<T>::name_to_desc_and_idx(name, this->var_descs);
    return make_contiguous_variable_view(var_desc, &data_array[idx]);
  };
private:
  // Note that we need the data to be in a single location for the lifetime of
  // this object. As long as none of the methods change the size of the vector,
  // we should be fine, but it may be safer to use a different type (that can
  // not dynamically reallocate) for internal array storage.
  std::vector<T> data_array;
};

}

#endif // SPAECIES_VARIABLE_ARRAY_HPP
