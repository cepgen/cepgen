#include "CepGen/Utils/Functional.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Core/Exception.h"

#include <muParser.h>

namespace cepgen
{
  namespace utils
  {
    class FunctionalMuParser : public Functional
    {
      public:
        explicit FunctionalMuParser( const ParametersList& );
        double eval( const std::vector<double>& ) const override;

      private:
        mu::Parser parser_;
    };

    FunctionalMuParser::FunctionalMuParser( const ParametersList& params ) :
      Functional( params )
    {
      try {
        for ( size_t i = 0; i < vars_.size(); ++i )
          parser_.DefineVar( vars_[i], &values_[i] );
        parser_.SetExpr( expression_ );
      } catch ( const mu::Parser::exception_type& e ) {
        throw CG_ERROR( "FunctionalMuParser" )
          << "Failed to define the function\n\t"
          << expression_ << "\n\t"
          << std::string( e.GetPos(), '-' )+"^" << "\n\t"
          << e.GetMsg();
      }
    }

    double
    FunctionalMuParser::eval( const std::vector<double>& x ) const
    {
      try {
        return parser_.Eval();
      } catch ( const mu::Parser::exception_type& e ) {
        throw CG_WARNING( "FunctionalMuParser" )
          << "Failed to evaluate the function\n\t"
          << expression_ << "\n\t"
          << std::string( e.GetPos(), '-' )+"^" << "\n\t"
          << e.GetMsg();
      }
    }
  }
}

REGISTER_FUNCTIONAL( "MuParser", FunctionalMuParser )
