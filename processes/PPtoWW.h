#ifndef PPtoWW_h
#define PPtoWW_h

#include "processes/GenericKTProcess.h"

/// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process using \f$k_T\f$-factorization approach
class PPtoWW : public GenericKTProcess
{
 public:
  PPtoWW();
  inline ~PPtoWW() {;}
  
 private:
  void PrepareKTKinematics();
  double ComputeJacobian();
  double ComputeKTFactorisedMatrixElement();
  void FillCentralParticlesKinematics();
  
  /// Minimal rapidity of the first outgoing W boson
  double fYmin;
  /// Maximal rapidity of the first outgoing W boson
  double fYmax;
  /// Rapidity of the first outgoing W boson
  double fY1;
  /// Rapidity of the first outgoing W boson
  double fY2;
  /// Transverse momentum difference for the two outgoing W bosons
  double _ptdiff;
  /// Azimuthal angle difference for the two outgoing W bosons
  double _phiptdiff;
  
  // first outgoing W boson
  Particle::Momentum fPw1;
  // second outgoing W boson
  Particle::Momentum fPw2;
};

#endif

