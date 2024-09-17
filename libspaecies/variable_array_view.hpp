#ifndef SPAECIES_VARIABLE_ARRAY_VIEW_HPP
#define SPAECIES_VARIABLE_ARRAY_VIEW_HPP

#include <cstring>
#include <initializer_list>
#include <memory>
#include <tuple>
#include <type_traits>

#include "contiguous_variable_view.hpp"
#include "exceptions.hpp"

namespace spaecies {

template <class T>
class VariableArrayView {
public:
  using NonConstT = std::remove_const_t<T>;
  VariableArrayView(std::initializer_list<VarDescPtr> var_descs)
    : VariableArrayView(std::vector<VarDescPtr>(var_descs)) {
  }
  VariableArrayView(const std::vector<VarDescPtr>& var_descs)
    : var_desc_vec(var_descs),
      data_size(calc_size(var_descs)),
      owning_ptr(new T[data_size]), // After C++20, can use make_shared here.
      data_ptr(owning_ptr.get()) {
  };
  VariableArrayView(const std::vector<VarDescPtr>& var_descs, T* data_ptr)
    : var_desc_vec(var_descs),
      data_size(calc_size(var_descs)),
      data_ptr(data_ptr) {
  };
  // Cast to const view.
  operator VariableArrayView<const T>() const {
    return VariableArrayView<const T>(var_desc_vec, data_size, owning_ptr, data_ptr);
  }
  // Required for access to private constructor used above.
  friend class VariableArrayView<NonConstT>;
  ContiguousVariableView<T> get_variable(const std::string& name) const {
    auto [var_desc, idx] = name_to_desc_and_idx(name, var_descs());
    return ContiguousVariableView(var_desc, data_ptr + idx);
  };
  inline T* data() const {
    return data_ptr;
  };
  inline const std::vector<VarDescPtr>& var_descs() const {
    return var_desc_vec;
  }
  inline std::size_t size() const {
    return data_size;
  }
  VariableArrayView<NonConstT> deep_copy() const {
    VariableArrayView<NonConstT> copy(var_descs());
    copy_data_to_location(copy.data());
    return copy;
  }
protected:
  std::vector<VarDescPtr> var_desc_vec;
  std::size_t data_size;
  std::shared_ptr<T[]> owning_ptr;
  T* data_ptr;
  inline void copy_data_to_location(NonConstT* dest) const {
    std::memcpy(dest, data_ptr, data_size * sizeof(T));
  }
  static std::tuple<VarDescPtr, std::size_t> name_to_desc_and_idx(const std::string& name,
                                                                  const std::vector<VarDescPtr> var_descs) {
    std::size_t idx = 0;
    VarDescPtr var_desc_out;
    for (VarDescPtr var_desc : var_descs) {
      if (var_desc->name == name) {
        var_desc_out = var_desc;
        return {var_desc_out, idx};
      }
      idx += var_desc->size();
    }
    throw(VariableNotFoundException(name, "variable not found in variable array"));
  }
  static std::size_t calc_size(const std::vector<VarDescPtr> var_descs) {
    std::size_t my_size = 0;
    for (VarDescPtr var_desc : var_descs) {
      my_size += var_desc->size();
    }
    return my_size;
  }
private:
  // Raw constructor for casting purposes.
  VariableArrayView(const std::vector<VarDescPtr>& var_descs, std::size_t data_size,
                    std::shared_ptr<T[]> owning_ptr, T* data_ptr)
    : var_desc_vec(var_descs), data_size(data_size), owning_ptr(owning_ptr), data_ptr(data_ptr) {
  }
};

static_assert(std::is_copy_constructible_v<VariableArrayView<double>> && std::is_copy_assignable_v<VariableArrayView<double>>);
static_assert(std::is_move_constructible_v<VariableArrayView<double>> && std::is_move_assignable_v<VariableArrayView<double>>);

}

#endif // SPAECIES_VARIABLE_ARRAY_VIEW_HPP
