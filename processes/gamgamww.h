#ifndef _GAMGAMWW_H
#define _GAMGAMWW_H

#include "../include/process.h"

/**
 * @brief Computes the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process
 */
class GamGamWW : public Process
{
 public:
  GamGamWW();
  ~GamGamWW();
  double ComputeWeight();
 private:

};

#endif

