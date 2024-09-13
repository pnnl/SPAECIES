#ifndef SPAECIES_VARIABLE_VIEW_HPP
#define SPAECIES_VARIABLE_VIEW_HPP

#include "exceptions.hpp"
#include "variable_descriptor.hpp"

namespace spaecies {

template <class T>
class VariableView {
public:
  // Get the variable descriptor for this variable.
  inline const VarDescPtr& var_desc() const {
    return var_desc_ptr;
  }
  // Access a scalar value of the flattened variable.
  virtual T& operator[](std::size_t idx) = 0;
  virtual const T operator[](std::size_t idx) const = 0;
  // Get the variable size.
  std::size_t size() const {
    return var_desc()->size();
  }

protected:
  VarDescPtr var_desc_ptr;

  VariableView(const VarDescPtr& var_desc) : var_desc_ptr(var_desc) {
    VariableType expected_type = SPAECIES_TYPE<T>;
    if (var_desc_ptr->type != expected_type) {
      throw TypeMismatchException(spaecies_type_name(var_desc_ptr->type),
                                  spaecies_type_name(expected_type),
                                  "descriptor for " + var_desc_ptr->name
                                  + " does not match Variable type");
    }
  };
};

}

#endif // SPAECIES_VARIABLE_VIEW_HPP
