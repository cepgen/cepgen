#ifndef EPA_h
#define EPA_h

#include "Physics.h"

namespace EPA
{
  /// Photon generation mode
  enum PhotonMode {
    InvalidMode             = 0,
    WeizsackerWilliams      = 1,
    Transversal             = 2,
    TransversalLongitudinal = 3
  };
  extern Particle::Momentum fProton, fElectron;
  extern double fEPAmax, fYmin, fYmax;
  extern double fME2, fMP2, fS;
  /// 4-product of the electron/proton momenta
  extern double fElDotPr;
  /// Electron energy
  extern double fEEl;
  /// Mode of operation for the EPA
  extern PhotonMode fMode;
  extern PhysicsBoundaries fBoundaries;
  /// Define the incoming state and physics parameters before computation
  void InitialiseEPA( const Particle& el, const Particle& pr, const PhotonMode& mode, const PhysicsBoundaries& b );
  /// Prepare the limit values and constants before computation
  void PrepareEPA();
  /// Compute the outgoing electron and photon's kinematics
  /// \param[in] x1 First integration variable: y
  /// \param[in] x2 Second integration variable: Q2
  /// \param[in] x3 Third integration variable: theta/eta for the outgoing electron
  /// \param[out] out_ele Electron-like 4-momentum
  /// \param[out] out_gam Photon 4-momentum
  /// \param[out] lf Longitudinal fraction of the momentum
  bool EPA( double x1, double x2, double x3, double* q2, Particle::Momentum* out_ele, Particle::Momentum* out_gam, double* lf );
}

#endif
