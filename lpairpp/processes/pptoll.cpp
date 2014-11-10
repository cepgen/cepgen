#include "pptoll.h"

PPtoLL::PPtoLL()
{
  _name = "gamma,gamma->l+,l- (kT-factorization approach)";
}

PPtoLL::~PPtoLL()
{}

double
PPtoLL::ComputeWeight()
{
  return -1.;
}
