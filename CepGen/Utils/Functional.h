#ifndef CepGen_Utils_Functional_h
#define CepGen_Utils_Functional_h

#include "CepGen/Modules/NamedModule.h"

#include <vector>
#include <array>
#include <string>

namespace cepgen
{
  class ParametersList;
  namespace utils
  {
    /// A string-to-functional parser
    /// \author L. Forthomme <laurent.forthomme@cern.ch>
    /// \date 21 Aug 2017
    class Functional : public NamedModule<std::string>
    {
      public:
        /// Default constructor
        Functional( const ParametersList& params );
        /// Compute the functional for a given value of the variable (one-dimensional case)
        /// \param[in] x Variable value
        double operator()( double x ) const;
        /// Compute the functional for a given value of the variables
        /// \param[in] x Variables values
        double operator()( const std::vector<double>& x ) const;

        /// List of string variable names
        const std::vector<std::string>& variables() const { return vars_orig_; }
        /// String expression held by this functional parser
        const std::string& expression() const { return expression_orig_; }

      protected:
        /// Compute the functional for a given value of the variables
        /// \param[in] x Variables values
        virtual double eval( const std::vector<double>& x ) const = 0;

      private:
        std::vector<std::string> vars_orig_; ///< User-defined variables to be reached
        std::string expression_orig_; ///< User-defined expression

      protected:
        std::vector<std::string> vars_; ///< Computer-readable variable to be reached
        std::string expression_; ///< Computer-readable expression
        mutable std::vector<double> values_; ///< Last arguments list fed to the functional
    };
  }
}

#endif
