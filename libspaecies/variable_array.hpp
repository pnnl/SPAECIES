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
    return ContiguousVariable(var_desc_out, &data[idx]);
  };
  const std::size_t size;
private:
  const std::vector<VarDescPtr> var_descs;
  std::vector<T> data;
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
