#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ParticleProperties.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <iostream>

namespace CepGen
{
  namespace sf
  {
    const double Parameterisation::mp_ = part::mass( PDG::proton );
    const double Parameterisation::mp2_ = Parameterisation::mp_*Parameterisation::mp_;

    double
    Parameterisation::F1( double xbj, double q2 ) const
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
    Parameterisation::computeFL( double xbj, double q2, const sr::Parameterisation& ratio )
    {
      double r_error = 0.;
      computeFL( xbj, q2, ratio( xbj, q2, r_error ) );
    }

    void
    Parameterisation::computeFL( double xbj, double q2, double r )
    {
      const double tau = 4.*xbj*xbj*mp2_/q2;
      FL = F2 * ( 1.+tau ) * ( r/( 1.+r ) );
    }

    std::string
    Parameterisation::description() const
    {
      std::ostringstream os;
      os << type;
      return os.str();
    }

    std::ostream&
    operator<<( std::ostream& os, const Parameterisation& sf )
    {
      os << sf.description();
      if ( sf.old_vals_ != std::pair<double,double>() )
        os << " at (" << sf.old_vals_.first << ", " << sf.old_vals_.second << "): "
           << "F2 = " << sf.F2 << ", FL = " << sf.FL;
      return os;
    }
  }

  /// Human-readable format of a structure function type
  std::ostream&
  operator<<( std::ostream& os, const sf::Type& sf )
  {
    switch ( sf ) {
      case sf::Type::Invalid:             return os << "[INVALID]";
      case sf::Type::Electron:            return os << "electron";
      case sf::Type::ElasticProton:       return os << "elastic proton";
      case sf::Type::SuriYennie:          return os << "Suri-Yennie";
      case sf::Type::SzczurekUleshchenko: return os << "Szczurek-Uleshchenko";
      case sf::Type::FioreBrasse:         return os << "Fiore-Brasse";
      case sf::Type::ChristyBosted:       return os << "Christy-Bosted";
      case sf::Type::CLAS:                return os << "CLAS";
      case sf::Type::BlockDurandHa:       return os << "BDH";
      case sf::Type::ALLM91:              return os << "ALLM91";
      case sf::Type::ALLM97:              return os << "ALLM97";
      case sf::Type::GD07p:               return os << "GD07p";
      case sf::Type::GD11p:               return os << "GD11p";
      case sf::Type::Schaefer:            return os << "LUXlike";
      case sf::Type::MSTWgrid:            return os << "MSTW (grid)";
      case sf::Type::LHAPDF:              return os << "LHAPDF";
    }
    return os;
  }
}
