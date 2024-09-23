#ifndef RAINSHAFT_TYPES_HPP
#define RAINSHAFT_TYPES_HPP
#include <vector>

#include "spaecies.hpp"

// Commonly-used types for the rainshaft model, mainly intended to provide
// shortened names for SPAECIES types.
// State (prognostic variables)
using State = spaecies::State<double>;
// Constant state (i.e. state treated as input-only)
using StateConst = spaecies::State<const double>;
// Tendency
using Tendency = spaecies::Tendency<double>;
// Constant variable view (for input-only variables)
using VarConst = spaecies::ContiguousVariableView<const double>;
// List of variable descriptors
using VarDescList = std::vector<spaecies::VarDescPtr>;
// Mutable variable view (for "outputs")
using VarMut = spaecies::ContiguousVariableView<double>;

#endif // RAINSHAFT_TYPES_HPP
