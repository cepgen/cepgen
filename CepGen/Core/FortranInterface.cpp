#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#ifdef __cplusplus
extern "C" {
#endif
  /// Expose structure functions calculators to Fortran
  void
  cepgen_structure_functions_( int& sfmode, double& xbj, double& q2, double& f2, double& fl )
  {
    using namespace cepgen;
    static auto sf = strfun::Parameterisation::build( ParametersList().set<int>( "id", sfmode ) );
    const auto& val = ( *sf )( xbj, q2 );
    f2 = val.F2;
    fl = val.FL;
  }

  double
  cepgen_kt_flux_( int& fmode, double& x, double& kt2, int& sfmode, double& mx )
  {
    using namespace cepgen;
    static auto sf = strfun::Parameterisation::build( ParametersList().set<int>( "id", sfmode ) );
    return ktFlux( (KTFlux)fmode, x, kt2, *sf, mx );
  }

  double
  cepgen_kt_flux_hi_( int& fmode, double& x, double& kt2, int& a, int& z )
  {
    using namespace cepgen;
    return ktFlux( (KTFlux)fmode, x, kt2, HeavyIon{ (unsigned short)a, (Element)z } );
  }

  double
  cepgen_particle_mass_( int& pdg_id )
  {
    try {
      return cepgen::PDG::get()( (cepgen::pdgid_t)pdg_id ).mass;
    } catch ( const cepgen::Exception& e ) {
      e.dump();
      exit( 0 );
    }
  }

  double
  cepgen_particle_charge_( int& pdg_id )
  {
    try {
      return cepgen::PDG::get()( (cepgen::pdgid_t)pdg_id ).charge/3.;
    } catch ( const cepgen::Exception& e ) {
      e.dump();
      exit( 0 );
    }
  }
#ifdef __cplusplus
}
#endif
