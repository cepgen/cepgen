#ifndef CepGen_Physics_PhotonFluxes_h
#define CepGen_Physics_PhotonFluxes_h

#include "Particle.h"
#include "FormFactors.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  /// Incoming parton fluxes
  namespace Fluxes
  {
    /// List of fluxes for incoming photons
    namespace Photon
    {
      /// Get the elastic flux to be expected at a given x_bjorken / kT
      /// \param[in] x Bjorken x
      /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\mathrm{T}}^2\f$ of the incoming photon
      double ProtonElastic( double x, double kt2 );

      /// Get the inelastic flux to be expected at a given x_bjorken / kT
      /// \param[in] x Bjorken x
      /// \param[in] kt2 Transverse 2-momentum \f$\mathbf{q}_{\mathrm{T}}^2\f$ of the incoming photon
      /// \param[in] mx Outgoing diffractive proton mass
      double ProtonInelastic( double x, double kt2, double mx );
    }
  }
}

#endif
