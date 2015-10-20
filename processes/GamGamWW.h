#ifndef _GAMGAMWW_H
#define _GAMGAMWW_H

#include "GenericProcess.h"

/**
 * @brief Computes the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process
 */
class GamGamWW : public GenericProcess
{
 public:
  GamGamWW();
  ~GamGamWW();
  double ComputeWeight();
 private:

};

#endif

