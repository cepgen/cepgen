#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/FormFactors.h"

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
    static auto sf = strfun::StructureFunctionsHandler::get().build( sfmode );
    const auto& val = ( *sf )( xbj, q2 );
    f2 = val.F2;
    fl = val.FL;
  }

  /// Compute a \f$k_{\rm T}\f$-dependent flux for single nucleons
  /// \param[in] fmode Flux mode (see cepgen::KTFlux)
  /// \param[in] x Fractional momentum loss
  /// \param[in] kt2 The \f$k_{\rm T}\f$ transverse momentum norm
  /// \param[in] sfmode Structure functions set for dissociative emission
  /// \param[in] mx Diffractive state mass for dissociative emission
  double
  cepgen_kt_flux_( int& fmode, double& x, double& kt2, int& sfmode, double& mx )
  {
    using namespace cepgen;
    static auto ff = ff::FormFactorsHandler::get().build( ff::Model::StandardDipole ); // use another argument for the modelling?
    ff->setStructureFunctions( strfun::StructureFunctionsHandler::get().build( sfmode ) );
    return ktFlux( (KTFlux)fmode, x, kt2, *ff, mx );
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
#ifdef __cplusplus
}
#endif
