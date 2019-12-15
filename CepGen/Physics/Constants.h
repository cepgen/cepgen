#ifndef CepGen_Physics_Constants_h
#define CepGen_Physics_Constants_h

#include <math.h>

namespace cepgen
{
  /// List of physical constants useful that may be used for the matrix element definition
  namespace constants
  {
    /// Electromagnetic coupling constant \f$\alpha_{\rm em}=\frac{e^2}{4\pi\epsilon_0\hbar c}\f$
    constexpr double ALPHA_EM = 1./137.035999;
    /// Electromagnetic charge (~0.303 in natural units)
    constexpr double G_EM = sqrt( 4.*M_PI*ALPHA_EM );
    /// Strong coupling constant \f$\alpha_{\rm QCD}\f$
    constexpr double ALPHA_QCD = 0.1184; // at the Z pole
    /// Conversion factor between GeV\f$^{-2}\f$ and barn
    /// i.e. \f$\hbar^2 c^2\f$ in GeV\f$^{-2}\f$.
    constexpr double GEVM2_TO_PB = 0.389351824e9;
    constexpr double SCONSTB = 2.1868465e10; // 1.1868465e10;
  }
}

#endif
