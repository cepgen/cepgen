#ifndef Physics_PhotonFluxes_h
#define Physics_PhotonFluxes_h

#include "Particle.h"
#include "FormFactors.h"
#include "core/Exception.h"

extern "C"
{
  extern void grv95lo_( float&, float&, float&, float&, float&, float&, float&, float& ); // xpart,q2part,uv,dv,us,ds,ss,wg
}

namespace PhotonFluxes
{
  /// Get the elastic flux to be expected at a given x_bjorken / kT
  double ProtonElastic(double x_, double kt2_);

#ifdef GRVPDF
  /// Get the inelastic flux to be expected at a given x_bjorken / kT
  double ProtonInelastic(double x_, double kt2_, double mx_);
#else
  inline double ProtonInelastic(double x_, double kt2_, double mx_) {
    InError( "Inelastic flux cannot be computed as GRV PDF set is not linked to this instance!" );
    exit(0);
  }
#endif
}

#endif
