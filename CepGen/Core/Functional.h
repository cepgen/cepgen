#ifndef CepGen_Core_Functional_h
#define CepGen_Core_Functional_h

#include <vector>
#include <array>
#include "CepGen/Core/utils.h"
#include "CepGen/Core/Exception.h"

#if defined MUPARSER
# include <muParser.h>
#elif defined EXPRTK
# include <exprtk.hpp>
#endif

using std::string;

namespace cepgen
{
  namespace utils
  {
    /// \brief A string-to-functional parser
    /// \tparam N Number of arguments
    /// \author L. Forthomme <laurent.forthomme@cern.ch>
    /// \date 21 Aug 2017
    template<size_t N>
    class Functional
    {
      public:
        /// Default constructor
        Functional() = default;
        /// Copy constructor
        Functional( const Functional& rhs ) :
#ifdef EXPRTK
          symbols_( rhs.symbols_ ), expr_( rhs.expr_ ),
#endif
          values_( rhs.values_ ), vars_( rhs.vars_ ), expression_( rhs.expression_ ) {
          initialise();
        }
        /// \brief Build a parser from an expression and a variables list
        /// \param[in] expr Expression to parse
        /// \param[in] vars List of variables to parse
        Functional( const std::string& expr, const std::vector<std::string>& vars ) :
          vars_( vars ), expression_( expr ) {
          initialise();
        }
        /// \brief Compute the functional for a given value of the variable (N=1 case)
        /// \param[in] x Variable value
        double eval( double x ) const {
          static_assert( N == 1, "This function only works with single-dimensional functions" );
          return eval( std::array<double,1>{ x } );
        }
        /// \brief Compute the functional for a given value of the variables
        /// \param[in] x Variables values
        double eval( const std::array<double,N>& x ) const {
#if defined MUPARSER
          values_ = x;
          try {
            return parser_.Eval();
          } catch ( const mu::Parser::exception_type& e ) {
            throw CG_WARNING( "Functional" )
              << "Failed to evaluate the function\n\t"
              << expression_ << "\n\t"
              << std::string( e.GetPos(), '-' )+"^" << "\n\t"
              << e.GetMsg();
          }
#elif defined EXPRTK
          values_ = x;
          return expr_.value();
#else
          throw CG_FATAL( "Functional" )
            << "Neither exprtk nor muParser are linked to this program.\n\t"
            << "The math evaluator is hence disabled!";
#endif
        }

      private:
        void initialise() {
          if ( vars_.size() != values_.size() )
            throw CG_FATAL( "Functional" )
              << "Number of values should match exactly the number of variables!";
#if defined MUPARSER
          try {
            for ( unsigned short i = 0; i < vars_.size(); ++i )
              parser_.DefineVar( vars_[i], &values_[i] );
            parser_.SetExpr( expression_ );
          } catch ( const mu::Parser::exception_type& e ) {
            std::ostringstream os;
            for ( unsigned short i = 0; i < e.GetPos(); ++i )
              os << "-"; os << "^";
            throw CG_WARNING( "Functional" )
              << "Failed to define the function\n\t"
              << expression_ << "\n\t"
              << std::string( e.GetPos(), '-' )+"^" << "\n\t"
              << e.GetMsg();
          }
#elif defined EXPRTK
          for ( unsigned short i = 0; i < vars_.size(); ++i )
            symbols_.add_variable( vars_[i], values_[i] );
          symbols_.add_constants();
          expr_.register_symbol_table( symbols_ );
          parser_.compile( expression_, expr_ );
#else
          throw CG_FATAL( "Functional" )
            << "Neither exprtk nor muParser are linked to this program.\n\t"
            << "The math evaluator is hence disabled!";
#endif
        }
#if defined MUPARSER
        mu::Parser parser_;
#elif defined EXPRTK
        exprtk::symbol_table<double> symbols_;
        exprtk::expression<double> expr_;
        exprtk::parser<double> parser_;
#endif
        mutable std::array<double,N> values_;
        std::vector<std::string> vars_;
        std::string expression_;
    };
  }
}

#endif
