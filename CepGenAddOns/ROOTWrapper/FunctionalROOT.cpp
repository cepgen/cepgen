#include "CepGen/Utils/Functional.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Core/Exception.h"

#include "TFormula.h"

namespace cepgen
{
  namespace utils
  {
    class FunctionalROOT : public Functional
    {
      public:
        FunctionalROOT( const ParametersList& params ) :
          Functional( params ) {
          for ( size_t i = 0; i < vars_.size(); ++i )
            func_.AddVariable( vars_[i], 0. );
          if ( func_.Compile( expression_.c_str() ) != 0 )
            throw CG_ERROR( "FunctionalROOT" )
              << "Failed to define the function\n\t"
              << expression_;
        }
        double eval( const std::vector<double>& x ) const override {
          if ( !func_.IsValid() )
            throw CG_WARNING( "FunctionalROOT" )
              << "Cannot evaluate the invalid function at " << x << ".";
          return func_.EvalPar( values_.data() );
        }

      private:
        TFormula func_;
    };
  }
}

REGISTER_FUNCTIONAL( "ROOT", FunctionalROOT )
