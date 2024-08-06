#ifndef DERIVATIVES_HPP
#define DERIVATIVES_HPP

#include <array>
#include <tuple>

template <std::size_t I>
using Grad = std::array<double, I>;

template <std::size_t I>
using ValGrad = std::tuple<double, Grad<I>>;

#endif