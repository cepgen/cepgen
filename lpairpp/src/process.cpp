#include "process.h"

Process::Process() :
  _ndim(0), _point_set(false),
  _setin(false), _setout(false), _setkin(false)
{
  // This is where the particles will be stored
  _ev = new Event();
  _name = "<invalid process>";
}

Process::~Process()
{
  if (_point_set) {
    delete[] _x;
  }
  delete _ev;
}

void
Process::SetPoint(const unsigned int ndim_,double x_[])
{
  // Number of dimensions on which the integration will be performed
  _ndim = ndim_;

  // Phase space coordinate becomes a protected attribute
  _x = new double[ndim_];

  std::copy(x_, x_+ndim_, _x);
  
  std::cout << std::setprecision(16);
#ifdef DEBUG
  DumpPoint();
#endif
  _point_set = true;
}

void
Process::DumpPoint()
{
  std::cout << "[GamGam::DumpPoint] number of integration parameters : " << _ndim << std::endl;
  for (unsigned int i=0; i<(unsigned int)_ndim; i++) {
    std::cout << "  x[" << i << "] = " << _x[i] << std::endl;
  }
}
