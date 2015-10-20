#include "GamGamWW.h"

GamGamWW::GamGamWW() : GenericProcess("gamma,gamma->W+,W-")
{}

GamGamWW::~GamGamWW()
{}

double
GamGamWW::ComputeWeight()
{
  return -1.;
}
