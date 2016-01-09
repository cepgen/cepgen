#ifndef GamGamWW_h
#define GamGamWW_h

#include "../include/GenericProcess.h"

/// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process
class GamGamWW : public GenericProcess
{
 public:
  GamGamWW();
  ~GamGamWW();
  double ComputeWeight();
 private:

};

#endif

