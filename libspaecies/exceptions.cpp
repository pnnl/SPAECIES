#include "exceptions.hpp"

namespace spaecies {

std::string dimension_not_found_string(const std::string& name, const std::string& details) {
  return "no dimension named '" + name + "' was found: " + details;
}

DimensionNotFoundException::DimensionNotFoundException(const std::string& name,
                                                       const std::string& details)
  : runtime_error(dimension_not_found_string(name, details)) {
}

DimensionNotFoundException::DimensionNotFoundException(const DimensionNotFoundException& other) noexcept
  : runtime_error(other) {
}

std::string type_mismatch_string(const std::string& type1, const std::string& type2,
                                 const std::string& details) {
  return "encountered dynamic type " + type1 + " when " + type2 + " was expected: " + details;
}

TypeMismatchException::
TypeMismatchException(const std::string& type1, const std::string& type2,
                      const std::string& details)
                      : std::runtime_error(type_mismatch_string(type1, type2, details)) {
}

TypeMismatchException::TypeMismatchException(const TypeMismatchException& other) noexcept
  : runtime_error(other) {
}

UnreachableException::UnreachableException(const std::string& what_arg)
  : std::logic_error(what_arg) {
}

UnreachableException::UnreachableException(const char* what_arg)
  : std::logic_error(what_arg) {
}

UnreachableException::UnreachableException(const UnreachableException& other) noexcept
  : std::logic_error(other) {
}

std::string variable_not_found_string(const std::string& name, const std::string& details) {
  return "no variable named '" + name + "' was found: " + details;
}

VariableNotFoundException::VariableNotFoundException(const std::string& name,
                                                     const std::string& details)
  : runtime_error(variable_not_found_string(name, details)) {
}

VariableNotFoundException::VariableNotFoundException(const VariableNotFoundException& other) noexcept
  : runtime_error(other) {
}

}
