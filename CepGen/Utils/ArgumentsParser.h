#ifndef CepGen_test_ArgumentsParser_h
#define CepGen_test_ArgumentsParser_h

#include <vector>
#include <string>
#include <map>

namespace cepgen
{
  /// A generic command line arguments parser
  /// \author Laurent Forthomme
  /// \date Nov 2017
  class ArgumentsParser
  {
    public:
      /// A parameter parsed from user's input
      struct Parameter {

        //----- parameters constructors

        //--- string
        /// A string parameter
        Parameter( std::string name, std::string description = "", std::string* var = nullptr, std::string default_value = "" );

        //----- unsigned/signed integers
        /// An unsigned integer parameter
        Parameter( std::string name, std::string description, unsigned int* var = nullptr, unsigned int default_value = 0 );
        /// An integer parameter
        Parameter( std::string name, std::string description, int* var = nullptr, int default_value = 0 );
        /// A boolean parameter
        Parameter( std::string name, std::string description, bool* var = nullptr, bool default_value = false );

        //--- floats
        /// A double-precision floating point parameter
        Parameter( std::string name, std::string description, double* var = nullptr, double default_value = -999.999 );

        //--- complex formats
        /// A vector of strings parameter
        Parameter( std::string name, std::string description, std::vector<std::string>* var = nullptr, std::vector<std::string> default_value = {} );
        /// A vector of integer parameter
        Parameter( std::string name, std::string description, std::vector<int>* var = nullptr, std::vector<int> default_value = {} );
        /// A vector of floating point parameter
        Parameter( std::string name, std::string description, std::vector<double>* var = nullptr, std::vector<double> default_value = {} );

        Parameter& parse();
        friend std::ostream& operator<<( std::ostream&, const Parameter& );

        //----- parameters attributes

	std::vector<std::string> name; ///< Computer-readable name
        std::string description; ///< User-friendly parameter description
        std::string value; ///< Value (or default value)
        bool optional;

        //----- parameters containers

        /// Pointer to a string variable possibly handled by this parameter
        std::string* str_variable;
        /// Pointer to a double-precision floating point variable possibly handled by this parameter
        double* float_variable;
        /// Pointer to an integer variable possibly handled by this parameter
        int* int_variable;
        /// Pointer to an unsigned integer variable possibly handled by this parameter
        unsigned int* uint_variable;
        /// Pointer to a boolean variable possibly handled by this parameter
        bool* bool_variable;
        /// Pointer to a vector of string variables possibly handled by this parameter
        std::vector<std::string>* vec_str_variable;
        /// Pointer to a vector of integer variables possibly handled by this parameter
        std::vector<int>* vec_int_variable;
        /// Pointer to a vector of floating point variables possibly handled by this parameter
        std::vector<double>* vec_float_variable;
      };

    public:
      /// Arguments parser constructor
      /// \param[in] argc Number of argument words (key + value) to be parsed
      /// \param[in] argv[] List of argument words (key + value) to be parsed
      ArgumentsParser( int argc, char* argv[] );
      /// Add a parameter required for the parser
      template<typename... Args> ArgumentsParser& addArgument( Args&&... args ) {
        params_.emplace_back( std::forward<Args>( args )... );
        params_.rbegin()->optional = false;
        return *this;
      }
      /// Add a non-mandatory parameters that can be parsed
      template<typename... Args> ArgumentsParser& addOptionalArgument( Args&&... args ) {
        params_.emplace_back( std::forward<Args>( args )... );
        params_.rbegin()->optional = true;
        return *this;
      }
      /// Associate the command-line arguments to parameters
      ArgumentsParser& parse();
      /// Read required and optional parameters
      std::string operator[]( std::string name ) const;
      /// Dump the list of arguments into the terminal
      void dump() const;
      /// Show usage
      void print_help() const;
      /// Return usage message
      std::string help_message() const;
      /// Are extra configuration flags found in arguments list?
      const std::vector<std::string>& extra_config() const { return extra_config_; }

    private:
      /// A collection of parameters
      typedef std::vector<Parameter> ParametersCollection;

      std::string command_name_;
      const ParametersCollection help_str_;
      const ParametersCollection config_str_;
      bool help_req_;
      ParametersCollection params_;
      std::vector<std::pair<std::string,std::string> > args_;
      std::vector<std::string> extra_config_;
  };
}

#endif
