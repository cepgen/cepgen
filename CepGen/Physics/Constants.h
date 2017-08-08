#ifndef Constants_h
#define Constants_h

#include <math.h>

/// List of physical constants useful that may be used for the matrix element definition
class Constants
{
 public:
  /// Electromagnetic coupling constant \f$\alpha_{em}=\frac{e^2}{4\pi\epsilon_0\hbar c}\f$
  static double alphaEM;
  static double alphaQCD;
  /// \f$\frac{1}{(\hbar c)^2}~[\mathrm b^{-1}]\f$?
  static double muBarn;
  /// Conversion factor between GeV^2 and barn
  static double GeV2toBarn;
  static double sconstb;
  static double alphaReduced;
};

#endif
