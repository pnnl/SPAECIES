#ifndef SPAECIES_EXCEPTIONS_HPP
#define SPAECIES_EXCEPTIONS_HPP

#include <stdexcept>
#include <string>

namespace spaecies {

class DimensionNotFoundException : public std::runtime_error {
public:
  // Construct with a string describing the dimension not found.
  DimensionNotFoundException(const std::string& name, const std::string& details);
  // Copy constructor with noexcept.
  DimensionNotFoundException(const DimensionNotFoundException& other) noexcept;

};

class TypeMismatchException : public std::runtime_error {
public:
  // Construct with strings describing the dynamic type, then the expected
  // type, then any further details that need to be reported.
  TypeMismatchException(const std::string& dyn_type, const std::string& expected_type,
                        const std::string& details);
  // Copy constructor with noexcept.
  TypeMismatchException(const TypeMismatchException& other) noexcept;

};

class UnreachableException : public std::logic_error {
public:
  UnreachableException(const std::string& what_arg);
  UnreachableException(const char* what_arg);
  UnreachableException(const UnreachableException& other) noexcept;
};

}

#endif // SPAECIES_EXCEPTIONS_HPP
