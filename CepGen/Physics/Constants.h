#ifndef CepGen_Physics_Constants_h
#define CepGen_Physics_Constants_h

#include <math.h>

namespace CepGen
{
  /// List of physical constants useful that may be used for the matrix element definition
  namespace constants
  {
    /// Electromagnetic coupling constant \f$\alpha_{\rm em}=\frac{e^2}{4\pi\epsilon_0\hbar c}\f$
    constexpr double alphaEM = 1./137.035;
    /// Strong coupling constant \f$\alpha_{\rm QCD}\f$
    constexpr double alphaQCD = 0.1184; // at the Z pole
    /// Conversion factor between GeV\f${}^2\f$ and barn
    constexpr double GeV2toBarn = 0.389351824e9; // 1.e4*(197.3271**2);
    constexpr double sconstb = 2.1868465e10; // 1.1868465e10;
  }
}

#endif
