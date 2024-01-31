#ifndef SPAECIES_VARIABLE_ARRAY_VIEW_HPP
#define SPAECIES_VARIABLE_ARRAY_VIEW_HPP

#include <initializer_list>
#include <tuple>

#include "contiguous_variable_view.hpp"
#include "exceptions.hpp"

namespace spaecies {

template<class T> class VariableArrayView;

template<class T>
VariableArrayView<T> make_variable_array_view(const std::vector<VarDescPtr>& var_descs, T* data_ptr) {
  return VariableArrayView(var_descs, data_ptr);
}

template<class T>
const VariableArrayView<T> make_variable_array_view(const std::vector<VarDescPtr>& var_descs, const T* data_ptr) {
  return VariableArrayView(var_descs, const_cast<T*>(data_ptr));
}

template <class T>
class VariableArrayView {
public:
  const std::vector<VarDescPtr> var_descs;
  const std::size_t size;
  VariableArrayView(const std::vector<VarDescPtr>& var_descs, T* data_ptr)
    : var_descs(var_descs),
      size(calc_size(var_descs)),
      data_ptr(data_ptr) {
  };
  ContiguousVariableView<T> get_variable(const std::string& name) {
    auto [var_desc, idx] = name_to_desc_and_idx(name, var_descs);
    return make_contiguous_variable_view(var_desc, data_ptr + idx);
  };
  const ContiguousVariableView<T> get_variable(const std::string& name) const {
    auto [var_desc, idx] = name_to_desc_and_idx(name, var_descs);
    return make_contiguous_variable_view(var_desc, data_ptr + idx);
  };
  inline T* data() {
    return data_ptr;
  };
  inline const T* data() const {
    return data_ptr;
  };
protected:
  T* data_ptr;
  static std::tuple<VarDescPtr, std::size_t> name_to_desc_and_idx(const std::string& name,
                                                                  const std::vector<VarDescPtr> var_descs) {
    std::size_t idx = 0;
    VarDescPtr var_desc_out;
    for (VarDescPtr var_desc : var_descs) {
      if (var_desc->name == name) {
        var_desc_out = var_desc;
        break;
      }
      idx += var_desc->size();
    }
    if (var_desc_out == nullptr) {
      throw(VariableNotFoundException(name, "variable not found in variable array"));
    }
    return std::tuple(var_desc_out, idx);
  }
  static std::size_t calc_size(const std::vector<VarDescPtr> var_descs) {
    std::size_t my_size = 0;
    for (VarDescPtr var_desc : var_descs) {
      my_size += var_desc->size();
    }
    return my_size;
  }
};

}

#endif // SPAECIES_VARIABLE_ARRAY_VIEW_HPP
