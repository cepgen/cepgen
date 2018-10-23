#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ParticleProperties.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <iostream>

namespace cepgen
{
  namespace strfun
  {
    const double Parameterisation::mp_ = particleproperties::mass( PDG::proton );
    const double Parameterisation::mp2_ = Parameterisation::mp_*Parameterisation::mp_;

    Parameterisation::Parameterisation( double f2, double fl ) :
      type( Type::Invalid ), F2( f2 ), FL( fl ), old_vals_({ 0., 0. }),
      r_ratio_( new sigrat::E143 )
    {}

    Parameterisation::Parameterisation( const Parameterisation& sf ) :
      type( sf.type ), F2( sf.F2 ), FL( sf.FL ), params_( sf.params_ ), old_vals_( sf.old_vals_ ),
      r_ratio_( sf.r_ratio_ )
    {}

    Parameterisation::Parameterisation( const ParametersList& params ) :
      type( (Type)params.get<int>( "id" ) ), F2( 0. ), FL( 0. ), params_( params ), old_vals_({ 0., 0. }),
      r_ratio_( sigrat::Parameterisation::build(
        params.get<ParametersList>( "sigmaRatio", ParametersList().set<int>( "id", (int)sigrat::Type::E143 ) )
      ) )
    {}

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
    Parameterisation::computeFL( double xbj, double q2 )
    {
      if ( !r_ratio_ )
        throw CG_FATAL( "StructureFunctions:FL" )
          << "Failed to retrieve a R-ratio calculator!";
      double r_error = 0.;
      computeFL( xbj, q2, (*r_ratio_)( xbj, q2, r_error ) );
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
  operator<<( std::ostream& os, const strfun::Type& sf )
  {
    switch ( sf ) {
      case strfun::Type::Invalid:             return os << "[INVALID]";
      case strfun::Type::Electron:            return os << "electron";
      case strfun::Type::ElasticProton:       return os << "elastic proton";
      case strfun::Type::SuriYennie:          return os << "Suri-Yennie";
      case strfun::Type::SzczurekUleshchenko: return os << "Szczurek-Uleshchenko";
      case strfun::Type::FioreBrasse:         return os << "Fiore-Brasse";
      case strfun::Type::ChristyBosted:       return os << "Christy-Bosted";
      case strfun::Type::CLAS:                return os << "CLAS";
      case strfun::Type::BlockDurandHa:       return os << "BDH";
      case strfun::Type::ALLM91:              return os << "ALLM91";
      case strfun::Type::ALLM97:              return os << "ALLM97";
      case strfun::Type::GD07p:               return os << "GD07p";
      case strfun::Type::GD11p:               return os << "GD11p";
      case strfun::Type::Schaefer:            return os << "LUXlike";
      case strfun::Type::MSTWgrid:            return os << "MSTW (grid)";
      case strfun::Type::Partonic:            return os << "Partonic";
    }
    return os;
  }
}
