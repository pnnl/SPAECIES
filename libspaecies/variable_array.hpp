#ifndef SPAECIES_VARIABLE_ARRAY_HPP
#define SPAECIES_VARIABLE_ARRAY_HPP

#include <initializer_list>

#include "contiguous_variable.hpp"

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
  ContiguousVariable<T> get_variable(const std::string& name) {
    std::size_t idx = 0;
    VarDescPtr var_desc_out;
    for (VarDescPtr var_desc : var_descs) {
      if (var_desc->name == name) {
        var_desc_out = var_desc;
        break;
      }
      idx += var_desc->size();
    }
    // SPS: need to cover case where variable is not found.
    return make_contiguous_variable(var_desc_out, &data[idx]);
  };
  // SPS: need to unit test const version and consolidate with non-const version.
  const ContiguousVariable<T> get_variable(const std::string& name) const {
    std::size_t idx = 0;
    VarDescPtr var_desc_out;
    for (VarDescPtr var_desc : var_descs) {
      if (var_desc->name == name) {
        var_desc_out = var_desc;
        break;
      }
      idx += var_desc->size();
    }
    // SPS: need to cover case where variable is not found.
    return make_contiguous_variable(var_desc_out, &data[idx]);
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
};

}

#endif // SPAECIES_VARIABLE_ARRAY_HPP
