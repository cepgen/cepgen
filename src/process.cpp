#include "process.h"

Process::Process() :
  fX(0), fNumDimensions(0), fIsPointSet(false),
  _setin(false), _setout(false), _setkin(false)
{
  // This is where the particles will be stored
  fEvent = new Event();
  _name = "<invalid process>";
}

Process::~Process()
{
  if (fIsPointSet) delete[] fX;
  delete fEvent;
}

void
Process::SetPoint(const unsigned int ndim_,double x_[])
{
  // Number of dimensions on which the integration will be performed
  fNumDimensions = ndim_;

  // Phase space coordinate becomes a protected attribute
  if (!fX) fX = new double[ndim_];

  std::copy(x_, x_+ndim_, fX);  
  fIsPointSet = true;
  if (kLoggingLevel>=Debug) DumpPoint();
}

void
Process::DumpPoint()
{
  std::ostringstream os;
  for (unsigned int i=0; i<(unsigned int)fNumDimensions; i++) {
    os << Form("  x(%2d) = %8.6f\n\t", i, fX[i]);
  }
  PrintInfo(Form("Number of integration parameters: %d\n\t"
                 "%s", fNumDimensions, os.str().c_str()));
}
