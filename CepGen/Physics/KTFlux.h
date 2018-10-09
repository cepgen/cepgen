#ifndef CepGen_Physics_KTFlux_h
#define CepGen_Physics_KTFlux_h

#include "CepGen/Physics/PDG.h"
#include <ostream>

namespace CepGen
{
  class StructureFunctions;
  class HeavyIon;
  /// Collection of fundamental constants for kT fluxes definition
  struct KTFluxParameters
  {
    static const double kMinKTFlux; ///< Minimal value taken for a kT-factorised flux
    static const double kMP; ///< Proton mass, un GeV/c\f${}^2\f$
    static const double kMP2; ///< Squared proton mass
  };
  /// Type of incoming partons fluxes
  enum class KTFlux
  {
    invalid = -1,
    P_Photon_Elastic = 0,
    P_Photon_Inelastic = 1,
    P_Photon_Inelastic_Budnev = 11,
    P_Gluon_KMR = 20,
    HI_Photon_Elastic = 100
  };
  /// Human version of the flux name
  std::ostream& operator<<( std::ostream&, const KTFlux& );
  /// \brief Compute the flux for a given parton x/kT
  /// \param[in] type Flux modelling
  /// \param[in] x Parton momentum fraction
  /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\mathrm{T}}^2\f$ of the incoming parton
  /// \param[in] sf Structure functions evaluator
  /// \param[in] mx Outgoing diffractive proton mass
  double ktFlux( const KTFlux& type, double x, double kt2, StructureFunctions& sf, double mx = KTFluxParameters::kMP );
  /// \brief Compute the flux (from heavy ion) for a given parton x/kT
  /// \param[in] type Flux modelling
  /// \param[in] x Parton momentum fraction
  /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\mathrm{T}}^2\f$ of the incoming parton
  /// \param[in] hi Heavy ion properties
  double ktFlux( const KTFlux& type, double x, double kt2, const HeavyIon& hi );
}

#endif
