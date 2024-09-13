#ifndef SPAECIES_CONTIGUOUS_VARIABLE_VIEW_HPP
#define SPAECIES_CONTIGUOUS_VARIABLE_VIEW_HPP

#include <type_traits>

#include "variable_view.hpp"

namespace spaecies {

template<class T> class ContiguousVariableView;

template<class T>
class ContiguousVariableView : public VariableView<T> {
public:
  // Construct a ContiguousVariableView from a description and pointer.
  ContiguousVariableView(const VarDescPtr& var_desc, T* data_ptr)
    : VariableView<T>(var_desc), data_ptr(data_ptr) {
  }
  // Allow conversion to const.
  operator ContiguousVariableView<const T>() {
    return ContiguousVariableView<const T>(this->var_desc(), data_ptr);
  }
  // Access a scalar value of the flattened variable.
  inline T& operator[](std::size_t idx) {
    return data_ptr[idx];
  }
  inline const T operator[](std::size_t idx) const {
    return data_ptr[idx];
  }

private:
  T *data_ptr;

};

static_assert(std::is_copy_constructible_v<ContiguousVariableView<double>> && std::is_copy_assignable_v<ContiguousVariableView<double>>);
static_assert(std::is_move_constructible_v<ContiguousVariableView<double>> && std::is_move_assignable_v<ContiguousVariableView<double>>);

}

#endif // SPAECIES_CONTIGUOUS_VARIABLE_VIEW_HPP
