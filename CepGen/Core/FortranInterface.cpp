#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Processes/FortranKTProcess.h"

#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"

#ifdef __cplusplus
extern "C" {
#endif
  /// Expose structure functions calculators to Fortran
  void
  cepgen_structure_functions_( int& sfmode, double& xbj, double& q2, double& f2, double& fl )
  {
    using namespace cepgen;
    static auto sf = strfun::StructureFunctionsFactory::get().build( sfmode );
    const auto& val = ( *sf )( xbj, q2 );
    f2 = val.F2;
    fl = val.FL;
  }

  /// Compute a \f$k_{\rm T}\f$-dependent flux for single nucleons
  /// \param[in] fmode Flux mode (see cepgen::KTFlux)
  /// \param[in] x Fractional momentum loss
  /// \param[in] kt2 The \f$k_{\rm T}\f$ transverse momentum norm
  /// \param[in] sfmode Structure functions set for dissociative emission
  /// \param[in] min Incoming particle mass
  /// \param[in] mout Diffractive state mass for dissociative emission
  double
  cepgen_kt_flux_( int& fmode, double& x, double& kt2, int& sfmode, double& min, double& mout )
  {
    using namespace cepgen;
    static auto sf = strfun::StructureFunctionsFactory::get().build( sfmode );
    return ktFlux( (KTFlux)fmode, x, kt2, *sf, min*min, mout*mout );
  }

  /// Compute a \f$k_{\rm T}\f$-dependent flux for heavy ions
  /// \param[in] fmode Flux mode (see cepgen::KTFlux)
  /// \param[in] x Fractional momentum loss
  /// \param[in] kt2 The \f$k_{\rm T}\f$ transverse momentum norm
  /// \param[in] a Mass number for the heavy ion
  /// \param[in] z Atomic number for the heavy ion
  double
  cepgen_kt_flux_hi_( int& fmode, double& x, double& kt2, int& a, int& z )
  {
    using namespace cepgen;
    return ktFlux( (KTFlux)fmode, x, kt2, HeavyIon{ (unsigned short)a, (Element)z } );
  }

  /// Mass of a particle, in GeV/c^2
  double
  cepgen_particle_mass_( int& pdg_id )
  {
    try {
      return cepgen::PDG::get().mass( (cepgen::pdgid_t)pdg_id );
    } catch ( const cepgen::Exception& e ) {
      e.dump();
      exit( 0 );
    }
  }

  /// Charge of a particle, in e
  double
  cepgen_particle_charge_( int& pdg_id )
  {
    try {
      return cepgen::PDG::get().charge( (cepgen::pdgid_t)pdg_id );
    } catch ( const cepgen::Exception& e ) {
      e.dump();
      exit( 0 );
    }
  }

  /// Colour factor of a particle
  double
  cepgen_particle_colour_( int& pdg_id )
  {
    try {
      return cepgen::PDG::get().colours( (cepgen::pdgid_t)pdg_id );
    } catch ( const cepgen::Exception& e ) {
      e.dump();
      exit( 0 );
    }
  }

  void
  cepgen_init_()
  {
    cepgen::Generator gen;
  }

  void
  cepgen_debug_( char* str, int size )
  {
    CG_DEBUG( "fortran_process" ) << std::string( str, size );
  }

  void
  cepgen_warning_( char* str, int size )
  {
    CG_WARNING( "fortran_process" ) << std::string( str, size );
  }

  void
  cepgen_error_( char* str, int size )
  {
    CG_ERROR( "fortran_process" ) << std::string( str, size );
  }

  void
  cepgen_fatal_( char* str, int size )
  {
    throw CG_FATAL( "fortran_process" ) << std::string( str, size );
  }

#ifdef __cplusplus
}
#endif
