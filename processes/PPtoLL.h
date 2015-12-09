#ifndef PPtoLL_h
#define PPtoLL_h

#include "../include/GenericProcess.h"

/**
 * @brief Computes the matrix element for a CE \f$\gamma\gamma\rightarrow \ell^+\ell^-\f$ process using \f$k_T\f$-factorization approach
 */
class PPtoLL : public GenericProcess
{
 public:
  PPtoLL();
  ~PPtoLL();
  int GetNdim(ProcessMode) const;
  double ComputeWeight();
  void FillKinematics(bool symmetrise_=false);
  void AddEventContent();
 private:
  double INCqqbar();
  
  double _q1t, _q2t, _phiq1t, _phiq2t;
  /**
   * Rapidity of the first outgoing lepton
   */
  double _y1;
  /**
   * Rapidity of the first outgoing lepton
   */
  double _y2;
  double _ptdiff, _phiptdiff;
  /**
   * Invariant mass of the first outgoing proton (or remnant), in GeV
   */
  double _mx;
  /**
   * Invariant mass of the second outgoing proton (or remnant), in GeV
   */
  double _my;
  
  double _px_0, _px_x, _px_y, _px_z;
  double _py_0, _py_x, _py_y, _py_z;
  double _pl1_0, _pl1_x, _pl1_y, _pl1_z;
  double _pl2_0, _pl2_x, _pl2_y, _pl2_z;
};

#endif

