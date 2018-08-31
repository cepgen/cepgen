#ifndef CepGen_Physics_KTFlux_h
#define CepGen_Physics_KTFlux_h

#include "CepGen/Physics/PDG.h"
#include <ostream>

namespace CepGen
{
  class StructureFunctions;
  struct KTFluxParameters
  {
    static const double kMinKTFlux, kMP, kMP2;
  };
  /// Type of incoming partons fluxes
  enum class KTFlux
  {
    P_Photon_Elastic = 0,
    P_Photon_Inelastic = 1,
    P_Photon_Inelastic_Budnev = 11,
  };
  std::ostream& operator<<( std::ostream&, const KTFlux& );
  /// Get the flux at a given parton x/kT
  /// \param[in] x Parton momentum fraction
  /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\mathrm{T}}^2\f$ of the incoming parton
  /// \param[in] sf Structure functions evaluator
  /// \param[in] mx Outgoing diffractive proton mass
  double ktFlux( const KTFlux& type, double x, double kt2, StructureFunctions& sf, double mx = KTFluxParameters::kMP );
}

#endif
