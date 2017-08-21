#ifndef CepGen_Cards_FunctionBuilder_h
#define CepGen_Cards_FunctionBuilder_h

#include <vector>
#include <array>
#include "CepGen/Core/utils.h"

#ifdef MATHEX
#include <mathex.h>
#endif

using std::string;

namespace CepGen
{
  /// A string-to-functional parser
  /// \author L. Forthomme <laurent.forthomme@cern.ch>
  /// \date 21 Aug 2017
  template<std::size_t N>
  class FunctionBuilder
  {
    public:
      /// Default constructor
      FunctionBuilder() {}
      /// Build a parser from an expression and a variables list
      /// \param[in] expr Expression to parse
      /// \param[in] vars List of variables to parse
      FunctionBuilder( const std::string& expr, const std::array<std::string,N>& vars )
#ifndef MATHEX
        InError( "MathEx is not linked to this program! the math evaluator is hence disabled!" );
#else
      : vars_( vars ) {
        parser_.expression( expr );
        for ( unsigned short i = 0; i < vars_.size(); ++i ) {
          parser_.addvar( vars_[i], &values_[i] );
        }
#endif
      }
      /// Compute the functional for a given value of the variable (N=1 case)
      /// \param[in] x Variable value
      double eval( double x ) {
#ifndef MATHEX
        return 1.0;
#else
        static_assert( N==1, "This function only works with single-dimensional functions" );
        values_[0] = x;
        double ret = 1.0;
        try { ret = parser_.eval(); } catch ( const smlib::mathex::error& e ) {
          throw Exception( __PRETTY_FUNCTION__, Form( "Failed to evaluate the function:\n\t%s", e.what() ), JustWarning );
        }
        return ret;
#endif
      }
      /// Compute the functional for a given value of the variables
      /// \param[in] x Variables values
      double eval( const std::array<double,N>& x ) {
        values_ = x;
        double ret = 1.0;
        try { ret = parser_.eval(); } catch ( const smlib::mathex::error& e ) {
          throw Exception( __PRETTY_FUNCTION__, Form( "Failed to evaluate the function:\n\t%s", e.what() ), JustWarning );
        }
        return ret;
      }
      /// Reference to the expression to be parsed
      std::string& expression() { return parser_.expression(); }

    private:
#ifdef MATHEX
      smlib::mathex parser_;
      std::array<std::string,N> vars_;
      std::array<double,N> values_;
#endif
  };
}

#endif
