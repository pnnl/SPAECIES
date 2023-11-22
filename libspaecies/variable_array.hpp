#ifndef SPAECIES_VARIABLE_ARRAY_HPP
#define SPAECIES_VARIABLE_ARRAY_HPP

#include <initializer_list>
#include <tuple>

#include "contiguous_variable_view.hpp"
#include "exceptions.hpp"

namespace spaecies {

template <class T>
class VariableArray {
public:
  VariableArray(std::initializer_list<VarDescPtr> var_descs)
    : VariableArray(std::vector<VarDescPtr>(var_descs)) {
  }
  VariableArray(const std::vector<VarDescPtr>& var_descs)
    : var_descs(var_descs),
      size(calc_size(var_descs)),
      data(size) {
  };
  ContiguousVariableView<T> get_variable(const std::string& name) {
    auto [var_desc, idx] = name_to_desc_and_idx(name, var_descs);
    return make_contiguous_variable_view(var_desc, &data[idx]);
  };
  const ContiguousVariableView<T> get_variable(const std::string& name) const {
    auto [var_desc, idx] = name_to_desc_and_idx(name, var_descs);
    return make_contiguous_variable_view(var_desc, &data[idx]);
  };
  const std::size_t size;
  const std::vector<VarDescPtr> var_descs;
  // SPS: need to make data a private, fixed-size array, and make some
  // kind of getter function to get the underlying pointer when really needed.
  std::vector<T> data;
private:
  static std::size_t calc_size(const std::vector<VarDescPtr> var_descs) {
    std::size_t my_size = 0;
    for (VarDescPtr var_desc : var_descs) {
      my_size += var_desc->size();
    }
    return my_size;
  }
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
};

}

#endif // SPAECIES_VARIABLE_ARRAY_HPP
