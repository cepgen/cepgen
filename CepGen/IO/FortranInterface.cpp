#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/Processes/GenericKTProcess.h"

#ifdef __cplusplus
extern "C" {
#endif

  void
  cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
  {
    using namespace CepGen;
    const StructureFunctions sf = StructureFunctionsBuilder::get( (StructureFunctions::Type)sfmode, q2, xbj );
    f2 = sf.F2;
    fl = sf.FL;
  }

  double
  cepgen_kt_flux_( int& fmode, double& kt2, double& x, int& sfmode, double& mx )
  {
    using namespace CepGen;
    using namespace CepGen::Process;
    return GenericKTProcess::flux( (GenericKTProcess::Flux)fmode, kt2, x,
                                   (StructureFunctions::Type)sfmode, mx );
  }

  double
  cepgen_kt_flux_hi_( int& fmode, double& kt2, double& x, int& a, int& z )
  {
    using namespace CepGen;
    using namespace CepGen::Process;
    return GenericKTProcess::flux( (GenericKTProcess::Flux)fmode, kt2, x,
                                   Kinematics::HeavyIon{ ( unsigned short )a, ( unsigned short )z } );
  }

  double
  cepgen_particle_mass_( int& pdg_id )
  {
    return CepGen::ParticleProperties::mass( (CepGen::PDG)pdg_id );
  }

  double
  cepgen_particle_charge_( int& pdg_id )
  {
    return CepGen::ParticleProperties::charge( pdg_id );
  }
#ifdef __cplusplus
}
#endif

