#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  namespace sigrat
  {
    std::shared_ptr<Parameterisation>
    Parameterisation::build( const ParametersList& params )
    {
      const Type& type = (Type)params.get<int>( "id" );
      switch ( type ) {
        case Type::E143:
          return std::make_shared<E143>( params );
        case Type::R1990:
          return std::make_shared<R1990>( params );
        case Type::CLAS:
          return std::make_shared<CLAS>( params );
        case Type::SibirtsevBlunden:
          return std::make_shared<SibirtsevBlunden>( params );
        default:
          throw CG_FATAL( "Rratio" ) << "Failed to build a R-ratio estimator with type=" << (int)type << "!";
      }
    }
  }

  /// Human-readable format of a R-ratio computation method
  std::ostream&
  operator<<( std::ostream& os, const sigrat::Type& sf )
  {
    switch ( sf ) {
      case sigrat::Type::E143:
        return os << "E143";
      case sigrat::Type::R1990:
        return os << "R1990";
      case sigrat::Type::CLAS:
        return os << "CLAS";
      case sigrat::Type::SibirtsevBlunden:
        return os << "SibirtsevBlunden";
    }
    return os;
  }
}

