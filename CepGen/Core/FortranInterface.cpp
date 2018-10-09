#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/ParticleProperties.h"

#include "CepGen/Core/Exception.h"

#ifdef __cplusplus
extern "C" {
#endif
  /// Expose structure functions calculators to Fortran
  void
  cepgen_structure_functions_( int& sfmode, double& xbj, double& q2, double& f2, double& fl )
  {
    using namespace CepGen;
    sf::Type sf_mode = (sf::Type)sfmode;

    CG_DEBUG( "cepgen_structure_functions" ) << sf_mode;

    static auto& val = sf::Parameterisation::build( sf_mode )->operator()( xbj, q2 );
    f2 = val.F2;
    fl = val.FL;
  }

  double
  cepgen_kt_flux_( int& fmode, double& x, double& kt2, int& sfmode, double& mx )
  {
    using namespace CepGen;
    static auto sf = sf::Parameterisation::build( (sf::Type)sfmode );
    return ktFlux(
      (KTFlux)fmode, x, kt2, *sf, mx );
  }

  double
  cepgen_kt_flux_hi_( int& fmode, double& x, double& kt2, int& a, int& z )
  {
    using namespace CepGen;
    return ktFlux(
      (KTFlux)fmode, x, kt2, HeavyIon{ (unsigned short)a, (Element)z } );
  }

  double
  cepgen_particle_mass_( int& pdg_id )
  {
    try {
      return CepGen::part::mass( (CepGen::PDG)pdg_id );
    } catch ( const CepGen::Exception& e ) {
      e.dump();
      exit( 0 );
    }
  }

  double
  cepgen_particle_charge_( int& pdg_id )
  {
    try {
      return CepGen::part::charge( pdg_id );
    } catch ( const CepGen::Exception& e ) {
      e.dump();
      exit( 0 );
    }
  }
#ifdef __cplusplus
}
#endif
