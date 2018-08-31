#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/ParticleProperties.h"

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

    static StructureFunctions& val = StructureFunctionsBuilder::get( sf_mode )->operator()( xbj, q2 );
    f2 = val.F2;
    fl = val.FL;
  }

  double
  cepgen_kt_flux_( int& fmode, double& x, double& kt2, int& sfmode, double& mx )
  {
    using namespace CepGen;
    static auto sf = StructureFunctionsBuilder::get( (SF::Type)sfmode );
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

namespace CepGen
{
  namespace Process
  {
    struct FortranKTProcessWrapper : GenericKTProcess
    {
      FortranKTProcessWrapper( const std::string& name, unsigned short parton_pdg, unsigned short outgoing_pdg ) :
        GenericKTProcess( name, std::string( "Fortran wrapped ")+name, 4,
                          { (ParticleCode)parton_pdg, (ParticleCode)parton_pdg },
                          { (ParticleCode)outgoing_pdg, (ParticleCode)outgoing_pdg } ) {}
      ~FortranKTProcessWrapper() {}
      double x[10];
    };
  }
}

extern "C"
{
  extern double test_weight_();
  extern struct {
  } test_params_;
}

#define add_kt_process( name, prefix, parton_pdg, outgoing_pdg ) \
  extern "C" {\
    extern double prefix##_weight_();\
    extern struct {\
    } prefix##_params_;\
  } \

add_kt_process( pptoll_fortran, pptoll, 22, 11 )
