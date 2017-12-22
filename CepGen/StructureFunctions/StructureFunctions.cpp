#include "StructureFunctions.h"
#include "CepGen/Physics/ParticleProperties.h"
#include <iostream>

namespace CepGen
{
  void
  StructureFunctions::computeFL( double xbj, double q2, double r )
  {
    const double mp2 = ParticleProperties::mass( Proton )*ParticleProperties::mass( Proton );
    const double tau = 4.*xbj*xbj*mp2/q2;
    FL = F2 * ( 1.+tau ) * ( r/( 1.+r ) );
  }

  /// Human-readable format of a structure function object
  std::ostream&
  operator<<( std::ostream& os, const StructureFunctions& sf )
  {
    return os << "F2 = " << sf.F2 << ", FL = " << sf.FL;
  }
  /// Human-readable format of a structure function type
  std::ostream&
  operator<<( std::ostream& os, const StructureFunctions::Type& sf )
  {
    switch ( sf ) {
      case StructureFunctions::Electron:            return os << "electron";
      case StructureFunctions::ElasticProton:       return os << "elastic proton";
      case StructureFunctions::SuriYennie:          return os << "Suri-Yennie";
      case StructureFunctions::SzczurekUleshchenko: return os << "Szczurek-Uleshchenko";
      case StructureFunctions::FioreBrasse:         return os << "Fiore-Brasse";
      case StructureFunctions::ChristyBosted:       return os << "Christy-Bosted";
      case StructureFunctions::CLAS:                return os << "CLAS";
      case StructureFunctions::BlockDurandHa:       return os << "BDH";
      case StructureFunctions::ALLM91:              return os << "ALLM;91";
      case StructureFunctions::ALLM97:              return os << "ALLM;97";
      case StructureFunctions::GD07p:               return os << "ALLM;GD07p";
      case StructureFunctions::GD11p:               return os << "ALLM;GD11p";
      case StructureFunctions::Schaefer:            return os << "Schaefer";
      case StructureFunctions::MSTWgrid:            return os << "MSTW (grid)";
    }
    return os;
  }
}
