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
  }
}
