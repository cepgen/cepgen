#ifndef CepGen_Physics_KTFlux_h
#define CepGen_Physics_KTFlux_h

#include "CepGen/Physics/PDG.h"
#include <iosfwd>

namespace cepgen
{
  namespace strfun { class Parameterisation; }
  class HeavyIon;
  /// Collection of fundamental constants for \f$k_{\rm T}\f$ fluxes definition
  struct KTFluxParameters
  {
    static const double kMinKTFlux; ///< Minimal value taken for a \f$\k_{\rm T}\f$-factorised flux
  };
  /// Type of incoming partons fluxes
  enum class KTFlux
  {
    invalid = -1, ///< Invalid flux
    P_Photon_Elastic = 0, ///< Elastic photon emission from proton
    P_Photon_Inelastic = 1, ///< Inelastic photon emission from proton
    P_Photon_Inelastic_Budnev = 11, ///< Inelastic photon emission from proton (Budnev flux approximation)
    P_Gluon_KMR = 20, ///< Inelastic gluon emission from proton (KMR flux modelling)
    HI_Photon_Elastic = 100 ///< Elastic photon emission from heavy ion (from Starlight \cite Klein:2016yzr)
  };
  /// Human version of the flux name
  std::ostream& operator<<( std::ostream&, const KTFlux& );
  /// \brief Compute the flux for a given parton \f$(x,k_{\rm T})\f$
  /// \param[in] type Flux modelling
  /// \param[in] x Parton momentum fraction
  /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\rm T}^2\f$ of the incoming parton
  /// \param[in] sf Structure functions evaluator
  /// \param[in] mi2 Incoming particle squared mass
  /// \param[in] mf2 Outgoing diffractive squared mass
  double ktFlux( const KTFlux& type, double x, double kt2, strfun::Parameterisation& sf, double mi2, double mf2 );
  /// \brief Compute the flux (from heavy ion) for a given parton \f$(x,k_{\rm T})\f$
  /// \param[in] type Flux modelling
  /// \param[in] x Parton momentum fraction
  /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\rm T}^2\f$ of the incoming parton
  /// \param[in] hi Heavy ion properties
  double ktFlux( const KTFlux& type, double x, double kt2, const HeavyIon& hi );
}

#endif
