#ifndef PPtoLL_h
#define PPtoLL_h

#include "core/GenericKTProcess.h"

/// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow \ell^+\ell^-\f$ process using \f$k_T\f$-factorization approach
class PPtoLL : public GenericKTProcess
{
 public:
  PPtoLL();
  inline ~PPtoLL() {;}
  
 private:
  void PrepareKTKinematics();
  double ComputeJacobian();
  /// \note IncQQbar in pptoll
  double ComputeKTFactorisedMatrixElement();
  void FillCentralParticlesKinematics();
  
  /// Minimal rapidity of the first outgoing lepton
  double fYmin;
  /// Maximal rapidity of the first outgoing lepton
  double fYmax;
  /// Rapidity of the first outgoing lepton
  double fY1;
  /// Rapidity of the first outgoing lepton
  double fY2;
  /// Transverse momentum difference for the two outgoing leptons
  double fPtDiff;
  /// Azimuthal angle difference for the two outgoing leptons
  double fPhiPtDiff;
  
  // first outgoing lepton
  Particle::Momentum fPl1;
  // second outgoing lepton
  Particle::Momentum fPl2;
};

#endif

