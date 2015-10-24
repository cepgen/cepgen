#include "include/MCGen.h"

int main()
{
  MCGen mg;
  Parameters *par = mg.parameters;

  par->in1p = par->in2p = 3500.;
  par->process = new PPtoLL;

  mg.parameters->Dump();
  double xsection, error_xsection;
  mg.ComputeXsection(&xsection, &error_xsection);

  return 0;
}
