#ifndef CepGen_Physics_TamingFunction_h
#define CepGen_Physics_TamingFunction_h

#include "CepGen/Core/Functional.h"
#include "CepGen/Core/utils.h"

namespace cepgen
{
  namespace utils
  {
    /// A collection of expression/variable with the associated functional
    struct TamingFunction
    {
      TamingFunction() {}
      TamingFunction( const std::string& var, const std::string& expr ) :
        var_orig( var ), expr_orig( expr ), var_safe( var ), expr_safe( expr ) {
        replace_all( var_safe, "(", "_" );
        replace_all( var_safe, ")", "_" );
        replace_all( expr_safe, var_orig, var_safe );
        function = Functional<1>( expr_safe, { { var_safe } } );
      }
      std::string var_orig, expr_orig, var_safe, expr_safe;
      Functional<1> function;
    };
    /// A collection of taming functions evaluator with helper classes
    class TamingFunctionsCollection : public std::vector<TamingFunction> {};
  }
}

#endif

