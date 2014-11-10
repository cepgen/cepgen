#ifndef _PPTOLL_H
#define _PPTOLL_H

#include "../include/process.h"

/**
 * @brief Computes the matrix element for a CE \f$\gamma\gamma\rightarrow \ell^+\ell^-\f$ process using \f$k_T\f$-factorization approach
 */
class PPtoLL : public Process
{
 public:
  PPtoLL();
  ~PPtoLL();
  double ComputeWeight();
 private:

};

#endif

