#ifndef CepGen_Core_Functional_h
#define CepGen_Core_Functional_h

#include <vector>
#include <array>
#include "CepGen/Core/utils.h"
#include "CepGen/Core/Exception.h"

#ifdef MUPARSER
#include <muParser.h>
#endif

using std::string;

namespace CepGen
{
  /// \brief A string-to-functional parser
  /// \tparam N Number of arguments
  /// \author L. Forthomme <laurent.forthomme@cern.ch>
  /// \date 21 Aug 2017
  template<std::size_t N>
  class Functional
  {
    public:
      /// Default constructor
      Functional() {}
      /// Copy constructor
      Functional( const Functional& rhs )
#ifdef MUPARSER
        : values_( rhs.values_ ), vars_( rhs.vars_ ), expression_( rhs.expression_ ) {
          for ( unsigned short i = 0; i < vars_.size(); ++i ) {
          parser_.DefineVar( vars_[i], &values_[i] );
        }
        parser_.SetExpr( expression_ );
      }
#else
      {}
#endif
      /// Build a parser from an expression and a variables list
      /// \param[in] expr Expression to parse
      /// \param[in] vars List of variables to parse
      Functional( const std::string& expr, const std::array<std::string,N>& vars ) : vars_( vars ), expression_( expr ) {
#ifdef MUPARSER
        for ( unsigned short i = 0; i < vars_.size(); ++i ) {
          parser_.DefineVar( vars_[i], &values_[i] );
        }
        parser_.SetExpr( expr );
#else
        InError( "muParser is not linked to this program! the math evaluator is hence disabled!" );
#endif
      }
      /// Compute the functional for a given value of the variable (N=1 case)
      /// \param[in] x Variable value
      double eval( double x ) const {
        static_assert( N==1, "This function only works with single-dimensional functions" );
        return eval( std::array<double,1>{ { x } } );
      }
      /// Compute the functional for a given value of the variables
      /// \param[in] x Variables values
      double eval( const std::array<double,N>& x ) const {
        double ret = 0.0;
#ifdef MUPARSER
        values_ = x;
        try { ret = parser_.Eval(); } catch ( const mu::Parser::exception_type& e ) {
          throw Exception( __PRETTY_FUNCTION__, Form( "Failed to evaluate the function:\n\t%s", e.GetMsg().c_str() ), JustWarning );
        }
#endif
        return ret;
      }
      /// Reference to the expression to be parsed
      //const std::string& expression() { return parser_.expression(); }

    private:
#ifdef MUPARSER
      mu::Parser parser_;
      mutable std::array<double,N> values_;
#endif
      std::array<std::string,N> vars_;
      std::string expression_;
  };
}

#endif
