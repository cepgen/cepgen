#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

namespace cepgen
{
  namespace utils
  {
    Functional::Functional( const ParametersList& params ) :
      NamedModule( params ),
      vars_orig_( params.get<std::vector<std::string> >( "variables" ) ),
      expression_orig_( params.get<std::string>( "expression" ) ),
      vars_( vars_orig_ ), expression_( expression_orig_ ),
      values_( vars_.size() )
    {
      for ( size_t i = 0; i < vars_.size(); ++i ) {
        replace_all( vars_.at( i ), "(", "_" );
        replace_all( vars_.at( i ), ")", "_" );
        replace_all( expression_, vars_orig_.at( i ), vars_.at( i ) );
      }
    }

    double
    Functional::operator()( double x ) const
    {
      if ( vars_orig_.size() != 1 )
        throw CG_FATAL( "Functional" )
          << "This function only works with single-dimensional functions!";
      return operator()( std::vector<double>{ x } );
    }

    double
    Functional::operator()( const std::vector<double>& x ) const
    {
      values_ = x;
      return eval( x );
    }
  }
}
