/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2024  Laurent Forthomme
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

#include <string>

#include "CepGen/Utils/Limits.h"

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
    inline ArgumentsParser& addArgument(Args&&... args) {
      params_.emplace_back(std::forward<Args>(args)...);
      params_.rbegin()->optional = false;
      return *this;
    }
    /// Add a non-mandatory parameters that can be parsed
    template <typename... Args>
    inline ArgumentsParser& addOptionalArgument(Args&&... args) {
      params_.emplace_back(std::forward<Args>(args)...);
      params_.rbegin()->optional = true;
      return *this;
    }
    ArgumentsParser& parse();                               ///< Associate command-line arguments to parameters
    std::string operator[](const std::string& name) const;  ///< Read required and optional parameters
    void dump() const;                                      ///< Dump the list of arguments into the terminal
    void print_help() const;                                ///< Show usage
    static void print_version();                            ///< Show version
    std::string help_message() const;                       ///< Usage message
    inline bool debugging() const { return debug_req_; }    ///< Is the debugging flag set?
    /// Extra configuration flags in arguments?
    inline const std::vector<std::string>& extra_config() const { return extra_config_; }

  private:
    /// A parameter parsed from user's input
    class Parameter {
    public:
      /// String parameter constructor
      Parameter(const std::string&, const std::string& = "", std::string* = nullptr, const std::string& = "");
      /// Unsigned integer parameter constructor
      Parameter(const std::string&, const std::string&, unsigned int* = nullptr, const unsigned int& = 0);
      /// Integer parameter constructor
      Parameter(const std::string&, const std::string&, int* = nullptr, const int& = 0);
      /// Boolean parameter constructor
      Parameter(const std::string&, const std::string&, bool* = nullptr, const bool& = false);
      /// Double-precision floating point parameter constructor
      Parameter(const std::string&, const std::string&, double* = nullptr, const double& = -999.999);
      /// Vector of strings parameter constructor
      Parameter(const std::string&,
                const std::string&,
                std::vector<std::string>* = nullptr,
                const std::vector<std::string>& = {});
      /// Vector of integer parameter constructor
      Parameter(const std::string&, const std::string&, std::vector<int>* = nullptr, const std::vector<int>& = {});
      /// Vector of floating point parameter constructor
      Parameter(const std::string&, const std::string&, std::vector<double>* = nullptr, const std::vector<double>& = {});
      /// Vector of floating point parameter constructor
      Parameter(const std::string&, const std::string&, Limits* = nullptr, const Limits& = Limits{});

      Parameter& parse();  ///< Cast the user input into a proper container value
      inline bool boolean() const { return bool_variable_ != nullptr; }  ///< Is the parameter a simple boolean?
      bool matches(const std::string&) const;  ///< Does the parameter name match a user-given argument?

      //----- parameters attributes

      std::vector<std::string> name;  ///< Computer-readable name
      std::string description;        ///< User-friendly parameter description
      std::string value;              ///< Value (or default value)
      bool optional{true};            ///< Flag to specify of the argument can be skipped from user input

    private:
      //----- parameters containers

      std::string* str_variable_{nullptr};                   ///< Pointer to a string variable
      double* float_variable_{nullptr};                      ///< Pointer to a double-precision floating point variable
      int* int_variable_{nullptr};                           ///< Pointer to an integer variable
      unsigned int* uint_variable_{nullptr};                 ///< Pointer to an unsigned integer variable
      bool* bool_variable_{nullptr};                         ///< Pointer to a boolean variable
      Limits* lim_variable_{nullptr};                        ///< Pointer to a limits object
      std::vector<std::string>* vec_str_variable_{nullptr};  ///< Pointer to a vector of string variables
      std::vector<int>* vec_int_variable_{nullptr};          ///< Pointer to a vector of integer variables
      std::vector<double>* vec_float_variable_{nullptr};     ///< Pointer to a vector of floating point variables
    };
    typedef std::vector<Parameter> ParametersCollection;  ///< A collection of parameters

    std::string command_name_;
    bool help_req_{false}, version_req_{false}, debug_req_{false};
    ParametersCollection params_;
    std::vector<std::pair<std::string, std::string> > args_;
    std::vector<std::string> add_ons_;
    std::vector<std::string> extra_config_;
  };
}  // namespace cepgen

#endif
