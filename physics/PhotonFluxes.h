#ifndef Physics_PhotonFluxes_h
#define Physics_PhotonFluxes_h

#include "Particle.h"
#include "FormFactors.h"
#include "core/Exception.h"

/// Get the elastic flux to be expected at a given x_bjorken / kT
double ElasticFlux(double x_, double kt2_);
#ifdef GRVPDF
/// Get the inelastic flux to be expected at a given x_bjorken / kT
double InelasticFlux(double x_, double kt2_, double mx_);
#else
inline double InelasticFlux(double x_, double kt2_, double mx_) {
  InError( "Inelastic flux cannot be computed as GRV PDF set is not linked to this instance!" );
  exit(0);
}
#endif

#endif
