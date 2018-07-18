#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ParticleProperties.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <iostream>

namespace CepGen
{
  const double StructureFunctions::mp_ = ParticleProperties::mass( PDG::Proton );
  const double StructureFunctions::mp2_ = StructureFunctions::mp_*StructureFunctions::mp_;

  double
  StructureFunctions::F1( double q2, double xbj ) const
  {
    if ( xbj == 0. || q2 == 0. ) {
      CG_ERROR( "StructureFunctions:F1" )
        << "Invalid range for Q² = " << q2 << " or xBj = " << xbj << ".";
      return 0.;
    }
    const double F1 = 0.5*( ( 1+4.*xbj*xbj*mp2_/q2 )*F2 - FL )/xbj;
    CG_DEBUG_LOOP( "StructureFunctions:F1" )
      << "F1 for Q² = " << q2 << ", xBj = " << xbj << ": " << F1 << "\n\t"
      << "(F2 = " << F2 << ", FL = " << FL << ").";
    return F1;
  }

  void
  StructureFunctions::computeFL( double q2, double xbj, const SF::SigmaRatio& ratio )
  {
    double r_error = 0.;
    computeFL( q2, xbj, ratio( q2, xbj, r_error ) );
  }

  void
  StructureFunctions::computeFL( double q2, double xbj, double r )
  {
    const double tau = 4.*xbj*xbj*mp2_/q2;
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
  operator<<( std::ostream& os, const SF::Type& sf )
  {
    switch ( sf ) {
      case SF::Type::Invalid:             return os << "[INVALID]";
      case SF::Type::Electron:            return os << "electron";
      case SF::Type::ElasticProton:       return os << "elastic proton";
      case SF::Type::SuriYennie:          return os << "Suri-Yennie";
      case SF::Type::SzczurekUleshchenko: return os << "Szczurek-Uleshchenko";
      case SF::Type::FioreBrasse:         return os << "Fiore-Brasse";
      case SF::Type::ChristyBosted:       return os << "Christy-Bosted";
      case SF::Type::CLAS:                return os << "CLAS";
      case SF::Type::BlockDurandHa:       return os << "BDH";
      case SF::Type::ALLM91:              return os << "ALLM;91";
      case SF::Type::ALLM97:              return os << "ALLM;97";
      case SF::Type::GD07p:               return os << "ALLM;GD07p";
      case SF::Type::GD11p:               return os << "ALLM;GD11p";
      case SF::Type::Schaefer:            return os << "Schaefer";
      case SF::Type::MSTWgrid:            return os << "MSTW (grid)";
      case SF::Type::GenericLHAPDF:       return os << "LHAPDF (generic)";
    }
    return os;
  }
}
