#ifndef SPAECIES_VARIABLE_ARRAY_VIEW_HPP
#define SPAECIES_VARIABLE_ARRAY_VIEW_HPP

#include <algorithm>
#include <cstring>
#include <initializer_list>
#include <memory>
#include <optional>
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
  VariableArrayView(std::vector<VarDescPtr> var_descs, T* data_ptr)
    : var_desc_vec(std::move(var_descs)),
      data_size(calc_size(var_desc_vec)),
      prognostic_data_size(calc_prognostic_size(var_desc_vec)),
      data_ptr(data_ptr) {
  };
  VariableArrayView(std::vector<VarDescPtr> var_descs)
    : var_desc_vec(std::move(var_descs)),
      data_size(calc_size(var_desc_vec)),
      prognostic_data_size(calc_prognostic_size(var_desc_vec)),
      owning_ptr(new T[data_size]), // After C++20, can use make_shared here.
      data_ptr(owning_ptr.get()) {
  };
  // Cast to const view.
  operator VariableArrayView<const T>() const {
    return VariableArrayView<const T>(var_desc_vec, data_size, prognostic_data_size, owning_ptr, data_ptr);
  }
  // Required for access to private constructor used above.
  friend class VariableArrayView<NonConstT>;
  std::optional<ContiguousVariableView<const T>> get_variable(const std::string& name) const {
    const std::optional<std::tuple<VarDescPtr, std::size_t>> desc_and_idx =
      name_to_desc_and_idx(name, var_descs());
    if (!desc_and_idx.has_value()) {
      return std::nullopt;
    }
    const auto [var_desc, idx] = desc_and_idx.value();
    return ContiguousVariableView<const T>{var_desc, data_ptr + idx};
  };
  std::optional<ContiguousVariableView<T>> get_variable(const std::string& name) {
    const std::optional<std::tuple<VarDescPtr, std::size_t>> desc_and_idx =
      name_to_desc_and_idx(name, var_descs());
    if (!desc_and_idx.has_value()) {
      return std::nullopt;
    }
    const auto [var_desc, idx] = desc_and_idx.value();
    return ContiguousVariableView<T>{var_desc, data_ptr + idx};
  };
  std::size_t get_idx(const std::string& name) const {
    return name_to_idx(name, var_descs());
  }
  inline T* data() {
    return data_ptr;
  };
  inline const T* data() const {
    return data_ptr;
  };
  inline const std::vector<VarDescPtr>& var_descs() const {
    return var_desc_vec;
  }
  inline std::size_t size() const {
    return data_size;
  }
  inline std::size_t prognostic_size() const {
    return prognostic_data_size;
  }
  VariableArrayView<NonConstT> deep_copy() const {
    VariableArrayView<NonConstT> copy(var_descs());
    copy_data_to_location(copy.data());
    return copy;
  }
protected:
  std::vector<VarDescPtr> var_desc_vec;
  std::size_t data_size;
  std::size_t prognostic_data_size;
  std::shared_ptr<T[]> owning_ptr;
  T* data_ptr;
  inline void copy_data_to_location(NonConstT* dest) const {
    std::copy_n(data_ptr, data_size, dest);
  }
  static std::size_t name_to_idx(const std::string& name, const std::vector<VarDescPtr>& var_descs) {
    std::size_t idx = 0;
    for (const VarDescPtr& var_desc : var_descs) {
      if (var_desc->name == name) {
        return idx;
      }
      idx += var_desc->size();
    }
    throw(VariableNotFoundException(name, "variable not found in variable array"));
  }
  static std::optional<std::tuple<VarDescPtr, std::size_t>> name_to_desc_and_idx(
      const std::string& name,
      const std::vector<VarDescPtr>& var_descs) {
    std::size_t idx = 0;
    for (const VarDescPtr& var_desc : var_descs) {
      if (var_desc->name == name) {
        return std::make_tuple(var_desc, idx);
      }
      idx += var_desc->size();
    }
    return std::nullopt;
  }
  static std::size_t calc_size(const std::vector<VarDescPtr> &var_descs) {
    std::size_t my_size = 0;
    for (const VarDescPtr& var_desc : var_descs) {
      my_size += var_desc->size();
    }
    return my_size;
  }
  static std::size_t calc_prognostic_size(const std::vector<VarDescPtr> &var_descs) {
    std::size_t my_size = 0;
    std::size_t i = 0;
    for (; i < var_descs.size() && var_descs[i]->usage == VariableUsage::Prognostic; i++) {
      my_size += var_descs[i]->size();
    }
    for (; i < var_descs.size(); i++) {
      if (var_descs[i]->usage == VariableUsage::Prognostic) {
        throw(InvalidArrayFormatException(var_descs[i]->name, "prognostic variables should be first"));
      }
    }

    return my_size;
  }
private:
  // Raw constructor for casting purposes.
  VariableArrayView(std::vector<VarDescPtr> var_descs, std::size_t data_size,
                    std::size_t prognostic_data_size, std::shared_ptr<T[]> owning_ptr, T* data_ptr)
    : var_desc_vec(std::move(var_descs)),
    data_size(data_size),
    prognostic_data_size(prognostic_data_size),
    owning_ptr(std::move(owning_ptr)),
    data_ptr(data_ptr) {
  }
};

static_assert(std::is_copy_constructible_v<VariableArrayView<double>> && std::is_copy_assignable_v<VariableArrayView<double>>);
static_assert(std::is_move_constructible_v<VariableArrayView<double>> && std::is_move_assignable_v<VariableArrayView<double>>);

}

#endif // SPAECIES_VARIABLE_ARRAY_VIEW_HPP
