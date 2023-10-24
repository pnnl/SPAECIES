#include "exceptions.hpp"

namespace spaecies {

std::string type_mismatch_string(const std::string& type1, const std::string& type2,
                                 const std::string& details) {
  return "encountered dynamic type " + type1 + " when " + type2 + " was expected: " + details;
}

TypeMismatchException::
TypeMismatchException(const std::string& type1, const std::string& type2,
                      const std::string& details)
                      : std::runtime_error(type_mismatch_string(type1, type2, details)) {
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

}
