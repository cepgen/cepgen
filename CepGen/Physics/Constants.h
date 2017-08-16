#ifndef CepGen_Physics_Constants_h
#define CepGen_Physics_Constants_h

#include <math.h>

namespace CepGen
{
  /// List of physical constants useful that may be used for the matrix element definition
  namespace Constants
  {
    /// Electromagnetic coupling constant \f$\alpha_\textrm{em}=\frac{e^2}{4\pi\epsilon_0\hbar c}\f$
    const double alphaEM = 1./137.035;
    /// Strong coupling constant \f$\alpha_\textrm{QCD}\f$
    const double alphaQCD = 0.1184; // at the Z pole
    /// Conversion factor between GeV^2 and barn
    const double GeV2toBarn = 3.89351824e8; // 1.e4*pow(197.3271, 2);
    const double sconstb = 2.1868465e10; // 1.1868465e10;
    const double alphaReduced = 0.5 * alphaEM / M_PI;
  }
}

#endif

