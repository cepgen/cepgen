#include "gamgamww.h"

GamGamWW::GamGamWW() : Process("gamma,gamma->W+,W-")
{}

GamGamWW::~GamGamWW()
{}

double
GamGamWW::ComputeWeight()
{
  return -1.;
}
