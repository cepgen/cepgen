/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Utils_ArgumentsParser_h
#define CepGen_Utils_ArgumentsParser_h

#include <map>
#include <string>
#include <vector>

namespace cepgen {
  /// A generic command line arguments parser
  /// \author Laurent Forthomme
  /// \date Nov 2017
  class ArgumentsParser {
  public:
    /// Arguments parser constructor
    /// \param[in] argc Number of argument words (key + value) to be parsed
    /// \param[in] argv[] List of argument words (key + value) to be parsed
    ArgumentsParser(int argc, char* argv[]);
    /// Add a parameter required for the parser
    template <typename... Args>
    ArgumentsParser& addArgument(Args&&... args) {
      params_.emplace_back(std::forward<Args>(args)...);
      params_.rbegin()->optional = false;
      return *this;
    }
    /// Add a non-mandatory parameters that can be parsed
    template <typename... Args>
    ArgumentsParser& addOptionalArgument(Args&&... args) {
      params_.emplace_back(std::forward<Args>(args)...);
      params_.rbegin()->optional = true;
      return *this;
    }
    /// Associate the command-line arguments to parameters
    ArgumentsParser& parse();
    /// Read required and optional parameters
    std::string operator[](std::string name) const;
    /// Dump the list of arguments into the terminal
    void dump() const;
    /// Show usage
    void print_help() const;
    /// Show version
    void print_version() const;
    /// Return usage message
    std::string help_message() const;
    /// Are extra configuration flags found in arguments list?
    const std::vector<std::string>& extra_config() const { return extra_config_; }

  private:
    /// A parameter parsed from user's input
    class Parameter {
    public:
      /// A string parameter constructor
      Parameter(std::string name,
                std::string description = "",
                std::string* var = nullptr,
                std::string default_value = "");
      /// An unsigned integer parameter constructor
      Parameter(std::string name, std::string description, unsigned int* var = nullptr, unsigned int default_value = 0);
      /// An integer parameter constructor
      Parameter(std::string name, std::string description, int* var = nullptr, int default_value = 0);
      /// A boolean parameter constructor
      Parameter(std::string name, std::string description, bool* var = nullptr, bool default_value = false);
      /// A double-precision floating point parameter constructor
      Parameter(std::string name, std::string description, double* var = nullptr, double default_value = -999.999);
      /// A vector of strings parameter constructor
      Parameter(std::string name,
                std::string description,
                std::vector<std::string>* var = nullptr,
                std::vector<std::string> default_value = {});
      /// A vector of integer parameter constructor
      Parameter(std::string name,
                std::string description,
                std::vector<int>* var = nullptr,
                std::vector<int> default_value = {});
      /// A vector of floating point parameter constructor
      Parameter(std::string name,
                std::string description,
                std::vector<double>* var = nullptr,
                std::vector<double> default_value = {});

      /// Cast the user input into a proper container value
      Parameter& parse();
      /// Is the parameter a simple boolean?
      inline bool boolean() const { return bool_variable_ != nullptr; }

      //----- parameters attributes

      std::vector<std::string> name;  ///< Computer-readable name
      std::string description;        ///< User-friendly parameter description
      std::string value;              ///< Value (or default value)
      bool optional;                  ///< Flag to specify of the argument can be skipped from user input

    private:
      //----- parameters containers

      /// Pointer to a string variable possibly handled by this parameter
      std::string* str_variable_;
      /// Pointer to a double-precision floating point variable possibly handled by this parameter
      double* float_variable_;
      /// Pointer to an integer variable possibly handled by this parameter
      int* int_variable_;
      /// Pointer to an unsigned integer variable possibly handled by this parameter
      unsigned int* uint_variable_;
      /// Pointer to a boolean variable possibly handled by this parameter
      bool* bool_variable_;
      /// Pointer to a vector of string variables possibly handled by this parameter
      std::vector<std::string>* vec_str_variable_;
      /// Pointer to a vector of integer variables possibly handled by this parameter
      std::vector<int>* vec_int_variable_;
      /// Pointer to a vector of floating point variables possibly handled by this parameter
      std::vector<double>* vec_float_variable_;
    };
    /// A collection of parameters
    typedef std::vector<Parameter> ParametersCollection;

    std::string command_name_;
    const ParametersCollection help_str_;
    const ParametersCollection version_str_;
    const ParametersCollection config_str_;
    bool help_req_{false}, version_req_{false};
    ParametersCollection params_;
    std::vector<std::pair<std::string, std::string> > args_;
    std::vector<std::string> extra_config_;
  };
}  // namespace cepgen

#endif
