#ifndef PPtoWW_h
#define PPtoWW_h

#include "../include/GenericKTProcess.h"

/// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process using \f$k_T\f$-factorization approach
class PPtoWW : public GenericKTProcess
{
 public:
  PPtoWW();
  ~PPtoWW();
  
  void PrepareKTKinematics();
  void FillKinematics(bool symmetrise_=false);
  
 private:
  double ComputeJacobian();
  /// Matrix element to be integrated
  /// \note IncQQbar in pptoll
  double ComputeKTFactorisedMatrixElement();
  
  double fYmin;
  double fYmax;
  /// Transverse virtuality of the first photon
  double _q1t;
  /// Transverse virtuality of the second photon
  double _q2t;
  /// Azimuthal rotation of the first photon transverse virtuality
  double _phiq1t;
  /// Azimuthal rotation of the first photon transverse virtuality
  double _phiq2t;
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

