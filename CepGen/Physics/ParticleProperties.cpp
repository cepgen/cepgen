#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  namespace particleproperties
  {
    double
    mass( const HeavyIon& hi )
    {
      if ( !hi )
        throw CG_FATAL( "mass" ) << "Invalid heavy ion: " << hi << "!";
      return (short)hi.Z*mass( PDG::proton );
    }

    double
    charge( const PDG& pdg_id )
    {
      switch ( pdg_id ) {
        case PDG::proton: case PDG::diffractiveProton:
          return +1.;
        case PDG::electron: case PDG::muon: case PDG::tau:
          return -1.;
        case PDG::down: case PDG::strange: case PDG::bottom:
          return -1./3;
        case PDG::up: case PDG::charm: case PDG::top:
          return +2./3;
        case PDG::W:
          return +1.;
        case PDG::piPlus: case PDG::KPlus: case PDG::DPlus:
          return +1.;
        default:
          return 0.;
      }
    }

    double
    charge( int id )
    {
      const short sign = id / abs( id );
      return sign * charge( (PDG)abs( id ) );
    }
  }
}
