#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/StructureFunctions/MSTWGrid.h"
#include "CepGen/Processes/GenericKTProcess.h"
#include "CepGen/Core/Exception.h"

#ifdef __cplusplus
extern "C" {
#endif
  void
  cepgen_structure_functions_( int& sfmode, double& xbj, double& q2, double& f2, double& fl )
  {
    using namespace CepGen;
    SF::Type sf_mode = (SF::Type)sfmode;

    CG_DEBUG( "cepgen_structure_functions" ) << sf_mode;

    static StructureFunctions& val = ( *StructureFunctionsBuilder::get( sf_mode ) )( xbj, q2 );
    f2 = val.F2;
    fl = val.FL;
  }

  double
  cepgen_kt_flux_( int& fmode, double& kt2, double& x, int& sfmode, double& mx )
  {
    using namespace CepGen;
    using namespace CepGen::Process;
    return GenericKTProcess::flux(
      (GenericKTProcess::Flux)fmode, kt2, x,
      *StructureFunctionsBuilder::get( (SF::Type)sfmode ), mx );
  }

  double
  cepgen_kt_flux_hi_( int& fmode, double& kt2, double& x, int& a, int& z )
  {
    using namespace CepGen;
    using namespace CepGen::Process;
    return GenericKTProcess::flux(
      (GenericKTProcess::Flux)fmode, kt2, x,
      HeavyIon{ (unsigned short)a, (unsigned short)z } );
  }

  double
  cepgen_particle_mass_( int& pdg_id )
  {
    try {
      return CepGen::ParticleProperties::mass( (CepGen::PDG)pdg_id );
    } catch ( const CepGen::Exception& e ) {
      e.dump();
      exit( 0 );
    }
  }

  double
  cepgen_particle_charge_( int& pdg_id )
  {
    try {
      return CepGen::ParticleProperties::charge( pdg_id );
    } catch ( const CepGen::Exception& e ) {
      e.dump();
      exit( 0 );
    }
  }
#ifdef __cplusplus
}
#endif

