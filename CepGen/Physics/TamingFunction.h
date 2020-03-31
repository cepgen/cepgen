#ifndef CepGen_Physics_TamingFunction_h
#define CepGen_Physics_TamingFunction_h

#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

namespace cepgen
{
  namespace utils
  {
    /// A pair of expression/variable with the associated functional
    struct TamingFunction
    {
      TamingFunction() {}
      /// Constructor for a taming function evaluator
      /// \param[in] var Variable to be tamed
      /// \param[in] expr String expression to define the taming
      TamingFunction( const std::string& var, const std::string& expr ) :
        var_orig( var ), expr_orig( expr ), var_safe( var ), expr_safe( expr ) {
        replace_all( var_safe, "(", "_" );
        replace_all( var_safe, ")", "_" );
        replace_all( expr_safe, var_orig, var_safe );
        function = Functional<1>( expr_safe, { { var_safe } } );
      }
      std::string var_orig; ///< User-defined variable to be tamed
      std::string expr_orig; ///< User-defined taming expression
      std::string var_safe; ///< Computer-readable variable to be tamed
      std::string expr_safe; ///< Computer-readable taming expression
      Functional<1> function; ///< Taming expression
    };
  }
}

#endif

