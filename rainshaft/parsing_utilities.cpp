#include "parsing_utilities.hpp"

namespace po = boost::program_options;

std::vector<std::string> split_string(const std::string &s, 
                                      const char delim) {
  std::vector<std::string> result;
  std::stringstream ss(s);
  std::string item;

  while (getline (ss, item, delim)) {
    result.push_back (item);
  }

  return result;
}

void parse_input_file(const std::string input_file) 
{
  std::ifstream settings_file(input_file);
  std::ofstream outfile("settings_filtered.ini");

  if (settings_file.is_open() && outfile.is_open()) {  
    std::string line; 

    while (std::getline(settings_file, line)) { 
      // Write each line to the output file
      if (line.find("=") != std::string::npos) {
        // Split the string
        std::vector<std::string> v = split_string(line, '=');

        // If variable on the line had a set value, i.e. key=value, write it to the filtered file
        if(v[1].find_first_not_of(' ') != std::string::npos) {
            outfile << line << "\n"; 
        } 
      } else {
        outfile << line << "\n"; 
      }
    }

    settings_file.close();  
    outfile.close();
  } else {
    throw std::logic_error("Unable to open settings file");
  }
}

void print_variables_map(po::variables_map vm) {
  for (const auto& it : vm) {
    std::cout << it.first.c_str() << " = ";
    auto& value = it.second.value();
    if (auto v = boost::any_cast<std::size_t>(&value))
      std::cout << *v;
    else if (auto v = boost::any_cast<int>(&value))
      std::cout << *v;
    else if (auto v = boost::any_cast<std::string>(&value))
      std::cout << *v;
    else if (auto v = boost::any_cast<double>(&value))
      std::cout << *v;
    else if (auto v = boost::any_cast<bool>(&value))
      std::cout << std::boolalpha << *v;
    else
      std::cout << "error printing variable value" << std::endl;
    
    if (it.second.defaulted())
      std::cout << " (defaulted)" << std::endl;
    else
      std::cout << std::endl;
  }
  std::cout << std::endl;
}