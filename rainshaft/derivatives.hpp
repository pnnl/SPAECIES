#ifndef DERIVATIVES_HPP
#define DERIVATIVES_HPP

#include <array>
#include <tuple>

template <std::size_t I>
using Grad = std::array<double, I>;

template <std::size_t I>
using ValGrad = std::tuple<double, Grad<I>>;

template <bool WithGrad, std::size_t I>
using Val = std::conditional_t<WithGrad, ValGrad<I>, double>;

constexpr double get_val(const double v) noexcept {
  return v;
}

template <std::size_t I>
constexpr double get_val(const ValGrad<I> &v) noexcept {
  return std::get<0>(v);
}

template <std::size_t I>
constexpr Grad<I> get_grad(const ValGrad<I> &v) noexcept {
  return std::get<1>(v);
}

#endif
