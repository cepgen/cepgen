#ifndef FormFactors_h
#define FormFactors_h

#include <math.h>

#include "physics/Constants.h"

/**
 * Compute the proton structure function (F.W Brasse et al., DESY 76/11 (1976),
 *   http://dx.doi.org/10.1016/0550-3213(76)90231-5)
 * @cite Brasse1976413
 */
bool PSF(double,double,double*,double*,double*);

/// Form factors collection (electric and magnetic parts)
struct FormFactors {
  /// Electric form factor
  double FE;
  /// Magnetic form factor
  double FM;
};
/// Trivial, spin-0 form factors (e.g. pion)
FormFactors TrivialFormFactors();

FormFactors ElasticFormFactors(double q2, double mi2);

FormFactors SuriYennieFormFactors(double q2, double mi2, double mf2);

/// Brasse et al. structure functions
/// \cite Brasse1976413
FormFactors FioreBrasseFormFactors(double q2, double mi2, double mf2);

#endif
