#ifndef PPtoLL_h
#define PPtoLL_h

#include "../include/GenericKTProcess.h"

/// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow \ell^+\ell^-\f$ process using \f$k_T\f$-factorization approach
class PPtoLL : public GenericKTProcess
{
 public:
  PPtoLL();
  ~PPtoLL();
  
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
  /// Rapidity of the first outgoing lepton
  double fY1;
  /// Rapidity of the first outgoing lepton
  double fY2;
  /// Transverse momentum difference for the two outgoing leptons
  double _ptdiff;
  /// Azimuthal angle difference for the two outgoing leptons
  double _phiptdiff;
  
  // first outgoing lepton
  Particle::Momentum fPl1;
  // second outgoing lepton
  Particle::Momentum fPl2;
};

#endif

