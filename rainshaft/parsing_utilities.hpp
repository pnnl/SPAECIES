#ifndef PARSING_UTILITIES_HPP
#define PARSING_UTILITIES_HPP

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>

namespace po = boost::program_options;

void parse_input_file(const std::string input_file);

std::vector<std::string> split_string(const std::string &s, const char delim);

void print_variables_map(const po::variables_map vm);

template <typename T>
void option_dependency(const po::variables_map &vm,
                       const std::string &for_what, 
                       const std::string &required_option, 
                       const T test_case={})
{
  if (test_case.empty()) {
    if (vm.count(for_what) && !vm[for_what].defaulted()) {
      if (vm.count(required_option) == 0 || vm[required_option].defaulted()) {
        throw std::logic_error(std::string("Option '") + for_what 
        + "'=" + test_case + " requires option '" + required_option + "'.");
      }
    }
  } else {
    if (vm.count(for_what) && !vm[for_what].defaulted()) {
      if (vm[for_what].as<T>() == test_case) {
        if (vm.count(required_option) == 0 || vm[required_option].defaulted()) {
          throw std::logic_error(std::string("Option '") + for_what 
          + "'=" + test_case + " requires option '" + required_option + "'.");
        }
      }
    }
  }
}

#endif // PARSING_UTILITIES_HPP