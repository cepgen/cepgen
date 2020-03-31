#ifndef CepGen_Core_Functional_h
#define CepGen_Core_Functional_h

#include <vector>
#include <array>
#include "CepGen/Utils/String.h"
#include "CepGen/Core/Exception.h"

#if defined FUNC_MUPARSER
# include <muParser.h>
#elif defined FUNC_EXPRTK
# include <exprtk.hpp>
#elif defined FUNC_ROOT
# include "TFormula.h"
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
#ifdef FUNC_EXPRTK
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
          values_ = x;
#if defined FUNC_MUPARSER
          try {
            return parser_.Eval();
          } catch ( const mu::Parser::exception_type& e ) {
            throw CG_WARNING( "Functional" )
              << "Failed to evaluate the function\n\t"
              << expression_ << "\n\t"
              << std::string( e.GetPos(), '-' )+"^" << "\n\t"
              << e.GetMsg();
          }
#elif defined FUNC_EXPRTK
          return expr_.value();
#elif defined FUNC_ROOT
          return func_.EvalPar( values_.data() );
#else
          throw CG_WARNING( "Functional" )
            << "Neither exprtk, muParser nor ROOT are linked to this program.\n\t"
            << "The formulas evaluator is hence disabled!";
#endif
        }
        /// String expression held by this functional parser
        const std::string& expression() const { return expression_; }

      private:
        void initialise() {
          if ( vars_.size() != values_.size() )
            throw CG_FATAL( "Functional" )
              << "Number of arguments (" << values_.size() << ") "
              << "should match exactly the number of variables\n"
              << "registered for this function " << "(" << vars_.size() << ")!";
#if defined FUNC_MUPARSER
          try {
            for ( size_t i = 0; i < vars_.size(); ++i )
              parser_.DefineVar( vars_[i], &values_[i] );
            parser_.SetExpr( expression_ );
          } catch ( const mu::Parser::exception_type& e ) {
            throw CG_WARNING( "Functional" )
              << "Failed to define the function\n\t"
              << expression_ << "\n\t"
              << std::string( e.GetPos(), '-' )+"^" << "\n\t"
              << e.GetMsg();
          }
#elif defined FUNC_EXPRTK
          for ( size_t i = 0; i < vars_.size(); ++i )
            symbols_.add_variable( vars_[i], values_[i] );
          symbols_.add_constants();
          expr_.register_symbol_table( symbols_ );
          parser_.compile( expression_, expr_ );
#elif defined FUNC_ROOT
          for ( size_t i = 0; i < vars_.size(); ++i )
            func_.AddVariable( vars_[i], 0. );
          if ( func_.Compile( expression_.c_str() ) != 0 )
            throw CG_WARNING( "Functional" )
              << "Failed to define the function\n\t"
              << expression_;
#else
          throw CG_WARNING( "Functional" )
            << "Neither exprtk, muParser nor ROOT are linked to this program.\n\t"
            << "The formulas evaluator is hence disabled!";
#endif
        }
#if defined FUNC_MUPARSER
        mu::Parser parser_;
#elif defined FUNC_EXPRTK
        exprtk::symbol_table<double> symbols_;
        exprtk::expression<double> expr_;
        exprtk::parser<double> parser_;
#elif defined FUNC_ROOT
        TFormula func_;
#endif
        mutable std::array<double,N> values_;
        std::vector<std::string> vars_;
        std::string expression_;
    };
  }
}

#endif
