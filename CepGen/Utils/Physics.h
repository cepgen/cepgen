#ifndef CepGen_Utils_Physics_h
#define CepGen_Utils_Physics_h

namespace cepgen
{
  namespace utils
  {
    /// Compute the diffractive mass from virtuality/Bjorken x
    double mX2( double xbj, double q2, double mp2 );
    /// Compute Bjorken x from virtuality/diffractive mass
    double xBj( double q2, double mp2, double mx2 );
  }
}

#endif
