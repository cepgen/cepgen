#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ParticleProperties.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <iostream>

namespace CepGen
{
  const double StructureFunctions::mp_ = ParticleProperties::mass( PDG::proton );
  const double StructureFunctions::mp2_ = StructureFunctions::mp_*StructureFunctions::mp_;

  double
  StructureFunctions::F1( double xbj, double q2 ) const
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
  StructureFunctions::computeFL( double xbj, double q2, const SF::SigmaRatio& ratio )
  {
    double r_error = 0.;
    computeFL( xbj, q2, ratio( xbj, q2, r_error ) );
  }

  void
  StructureFunctions::computeFL( double xbj, double q2, double r )
  {
    const double tau = 4.*xbj*xbj*mp2_/q2;
    FL = F2 * ( 1.+tau ) * ( r/( 1.+r ) );
  }

  std::string
  StructureFunctions::description() const
  {
    std::ostringstream os;
    os << type;
    return os.str();
  }

  /// Human-readable format of a structure function object
  std::ostream&
  operator<<( std::ostream& os, const StructureFunctions& sf )
  {
    os << sf.description();
    if ( sf.old_vals_ != std::pair<double,double>() )
      os << " at (" << sf.old_vals_.first << ", " << sf.old_vals_.second << "): "
         << "F2 = " << sf.F2 << ", FL = " << sf.FL;
    return os;
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
      case SF::Type::ALLM91:              return os << "ALLM91";
      case SF::Type::ALLM97:              return os << "ALLM97";
      case SF::Type::GD07p:               return os << "GD07p";
      case SF::Type::GD11p:               return os << "GD11p";
      case SF::Type::Schaefer:            return os << "LUXlike";
      case SF::Type::MSTWgrid:            return os << "MSTW (grid)";
      case SF::Type::LHAPDF:              return os << "LHAPDF";
    }
    return os;
  }
}
