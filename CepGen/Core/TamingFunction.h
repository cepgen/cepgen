#ifndef CepGen_Core_TamingFunction_h
#define CepGen_Core_TamingFunction_h

#include "Functional.h"
#include <unordered_map>

namespace CepGen
{
  /// A collection of expression/variable with the associated functional
  struct TamingFunction
  {
    TamingFunction( const std::string& var, const std::string& expr ) : variable( var ), expression( expr ), function( expr, { { variable } } ) {}
    std::string variable, expression;
    Functional<1> function;
  };
  /// A collection of taming functions evaluator with helper classes
  class TamingFunctionsCollection : public std::unordered_map<std::string, TamingFunction>
  {
    public:
      /// Insert a new variable/expression into the collection
      void add( const std::string& var, const std::string& expr ) { emplace( var, TamingFunction( var, expr ) ); }
      /// Does the collection handle a taming function for a given variable?
      bool has( const std::string& var ) const { return ( find( var ) != end() ); }
      /// Evaluate the taming function for a given variable at a given value
      double eval( const std::string& var, double x ) const {
        const auto& it = find( var );
        if ( it == end() ) return 1.0;
        return it->second.function.eval( x );
      }
      /// Dump a full list of taming functions handled
      void dump( std::ostream& os = *Logger::get().output ) const {
        os << "List of taming functions:\n";
        for ( const auto& it : *this )
          os << ">> \"" << it.second.expression << "\" applied on variable \"" << it.first << "\"\n";
      }
  };
}

#endif
